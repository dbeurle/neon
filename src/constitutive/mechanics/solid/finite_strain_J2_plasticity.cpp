
#include "finite_strain_J2_plasticity.hpp"

#include "exceptions.hpp"
#include "constitutive/internal_variables.hpp"
#include "numeric/float_compare.hpp"
#include "numeric/mechanics"

#include <range/v3/view/transform.hpp>
#include <range/v3/view/zip.hpp>

#include <unsupported/Eigen/MatrixFunctions>

#include <iostream>

namespace neon::mechanics::solid
{
finite_strain_J2_plasticity::finite_strain_J2_plasticity(std::shared_ptr<internal_variables_t>& variables,
                                                         json const& material_data)
    : small_strain_J2_plasticity(variables, material_data)
{
    variables->add(variable::second::hencky_strain_elastic);

    names.emplace("hencky_strain_elastic");

    // Add material tangent with the linear elasticity moduli
    variables->add(variable::fourth::tangent_operator,
                   consistent_tangent(1.0, matrix3::Zero(), matrix3::Zero(), C_e));
}

finite_strain_J2_plasticity::~finite_strain_J2_plasticity() = default;

void finite_strain_J2_plasticity::update_internal_variables(double)
{
    using namespace ranges;

    auto const shear_modulus = material.shear_modulus();

    // Extract the internal variables
    auto [deformation_gradients,
          log_strain_e_list,
          cauchy_stresses] = variables->get(variable::second::deformation_gradient,
                                            variable::second::hencky_strain_elastic,
                                            variable::second::cauchy_stress);

    auto const old_deformation_gradients = variables->get_old(variable::second::deformation_gradient);

    auto const J_list = variables->get(variable::scalar::DetF);

    // Retrieve the accumulated internal variables
    auto [accumulated_plastic_strains,
          von_mises_stresses] = variables->get(variable::scalar::effective_plastic_strain,
                                               variable::scalar::von_mises_stress);

    auto& tangent_operators = variables->get(variable::fourth::tangent_operator);

    auto const incremental_deformation_gradients = view::zip(deformation_gradients,
                                                             old_deformation_gradients)
                                                   | view::transform([](auto const& tpl) {
                                                         auto const& [F, F_old] = tpl;
                                                         return F * F_old.inverse();
                                                     });

    // Perform the update algorithm for each quadrature point
    for (std::size_t l{0}; l < deformation_gradients.size(); l++)
    {
        auto const& F_inc = incremental_deformation_gradients[l];
        auto const J = J_list[l];

        auto& cauchy_stress = cauchy_stresses[l];
        auto& accumulated_plastic_strain = accumulated_plastic_strains[l];
        auto& von_mises = von_mises_stresses[l];
        auto& log_strain_e = log_strain_e_list[l];

        // std::cout << "log strain elastic\n" << log_strain_e << std::endl;

        // Elastic trial deformation gradient
        matrix3 const B_e = (2.0 * log_strain_e).exp();

        // std::cout << "B elastic\n" << B_e << std::endl;

        // std::cout << "Finc\n" << F_inc << std::endl;

        // Elastic trial left Cauchy-Green deformation tensor
        matrix3 const B_e_trial = F_inc * B_e * F_inc.transpose();

        // std::cout << "B elastic trial\n" << B_e_trial << std::endl;

        // Trial Logarithmic elastic strain
        log_strain_e = 0.5 * B_e_trial.log();

        // std::cout << "log strain elastic\n" << log_strain_e << std::endl;

        // Elastic stress predictor
        cauchy_stress = compute_cauchy_stress(material.shear_modulus(), material.lambda(), log_strain_e)
                        / J;

        // Trial von Mises stress
        von_mises = von_mises_stress(cauchy_stress);

        // std::cout << "von Mises " << von_mises << std::endl;

        // Compute the initial estimate of the yield function for the material
        // and decide if the stress return needs to be computed
        if (auto const f = evaluate_yield_function(von_mises, accumulated_plastic_strain); f <= 0.0)
        {
            tangent_operators[l] = consistent_tangent(J, log_strain_e, cauchy_stress, C_e);
            continue;
        }

        auto const von_mises_trial = von_mises;

        // std::cout << "\nQUADRATURE POINT PLASTIC\n";

        // Compute the normal direction to the yield surface which remains
        // constant throughout the radial return method
        matrix3 const normal = deviatoric(cauchy_stress) / deviatoric(cauchy_stress).norm();

        // Initialise the plastic increment
        auto const plastic_increment = perform_radial_return(von_mises, accumulated_plastic_strain);

        // Plastic strain update
        // F_p = (plastic_increment * std::sqrt(3.0 / 2.0) * normal);

        // std::cout << ">>>Elastic strain before dec\n" << log_strain_e << std::endl;
        log_strain_e -= plastic_increment * std::sqrt(3.0 / 2.0) * normal;
        // std::cout << ">>>Elastic strain after dec\n" << log_strain_e << std::endl;

        cauchy_stress -= 2.0 * shear_modulus * plastic_increment * std::sqrt(3.0 / 2.0) * normal / J;

        von_mises = von_mises_stress(cauchy_stress);

        accumulated_plastic_strain += plastic_increment;

        // Use the elastoplastic tangent from infinitesimal strain theory
        matrix6 const D_ep = algorithmic_tangent(plastic_increment,
                                                 accumulated_plastic_strain,
                                                 von_mises_trial,
                                                 normal);

        // Compute the elastic-plastic tangent modulus for large strain
        tangent_operators[l] = consistent_tangent(J, log_strain_e, cauchy_stress, D_ep);
    }
}

matrix6 finite_strain_J2_plasticity::consistent_tangent(double const J,
                                                        matrix3 const& Be_trial,
                                                        matrix3 const& cauchy_stress,
                                                        matrix6 const& C)
{
    // Convert to Mandel notation so matrix multiplication == double dot operation
    matrix6 const D = voigt_to_mandel(C);
    matrix6 const L = voigt_to_mandel(compute_L(Be_trial));
    matrix6 const B = voigt_to_mandel(compute_B(Be_trial));

    matrix6 const H = compute_H(cauchy_stress);

    return 1.0 / (2.0 * J) * D * L * B - H;
}

std::tuple<vector3, std::array<matrix3, 3>, bool, std::array<int, 3>> finite_strain_J2_plasticity::
    compute_eigenvalues_eigenprojections(matrix3 const& X) const
{
    Eigen::SelfAdjointEigenSolver<matrix3> eigen_solver(X);

    if (eigen_solver.info() != Eigen::Success)
    {
        throw computational_error("Eigenvalue solver failed in finite plasticity routine\n");
    }

    vector3 const& x = eigen_solver.eigenvalues();
    matrix3 const& v = eigen_solver.eigenvectors();

    // Eigenprojections
    std::array<matrix3, 3> E = {{v.col(0) * v.col(0).transpose(),
                                 v.col(1) * v.col(1).transpose(),
                                 v.col(2) * v.col(2).transpose()}};

    bool is_repeated = false;
    std::array<int, 3> abc_ordering{{0, 1, 2}};

    // Need to check if there are repeated eigenvalues to machine precision
    // if (x.norm() < 1.0e-2 || (is_approx(x(0), x(1)) && is_approx(x(1), x(2))))
    if (is_approx(x(0), x(1)) && is_approx(x(1), x(2)))
    {
        std::cout << "We have repeated eigenvalues!\n" << x << std::endl;
        E[0] = matrix3::Identity();
        is_repeated = true;
        abc_ordering = {{-1, -1, -1}};
    }
    else if (!is_approx(x(0), x(1)) && is_approx(x(1), x(2)))
    {
        E[1] = matrix3::Identity() - E[0];
        is_repeated = true;
        abc_ordering = {{0, 1, 2}};
    }
    else if (!is_approx(x(2), x(0)) && is_approx(x(0), x(1)))
    {
        E[0] = matrix3::Identity() - E[2];
        is_repeated = true;
        abc_ordering = {{2, 0, 1}};
    }
    else if (!is_approx(x(1), x(2)) && is_approx(x(2), x(0)))
    {
        E[2] = matrix3::Identity() - E[1];
        is_repeated = true;
        abc_ordering = {{1, 2, 0}};
    }
    return {x, E, is_repeated, abc_ordering};
}

matrix6 finite_strain_J2_plasticity::derivative_tensor_log_unique(
    matrix3 const& Be_trial,
    vector3 const& x,
    std::array<matrix3, 3> const& E,
    std::array<int, 3> const& abc_ordering) const
{
    auto const [a, b, c] = abc_ordering;
    return std::log(x(a)) / ((x(a) - x(b)) * (x(a) - x(c)))
               * (dX2_dX(Be_trial) - (x(b) - x(c)) * Isym
                  - (2.0 * x(a) - x(b) - x(c)) * outer_product(E[a])
                  - (x(b) - x(c)) * (outer_product(E[b]) - outer_product(E[c])))
           + 1.0 / x(a) * outer_product(E[a]);
}

matrix6 finite_strain_J2_plasticity::compute_L(matrix3 const& Be_trial) const
{
    auto const [x, E, is_repeated, abc_ordering] = compute_eigenvalues_eigenprojections(Be_trial);

    if (!is_repeated)
    {
        return derivative_tensor_log_unique(Be_trial, x, E, {{0, 1, 2}})
               + derivative_tensor_log_unique(Be_trial, x, E, {{2, 0, 1}})
               + derivative_tensor_log_unique(Be_trial, x, E, {{1, 2, 0}});
    }
    else if (is_repeated && abc_ordering[0] >= 0)
    {
        // Derivative when there is one repeated eigenvalue
        auto const& [a, b, c] = abc_ordering;

        std::cout << "Eigenvalues\n" << x << std::endl;

        vector3 const y = x.array().log();

        std::cout << "log Eigenvalues\n" << y << std::endl;

        auto const s1 = (y(a) - y(c)) / std::pow(x(a) - x(c), 2) - 1.0 / (x(c) * (x(a) - x(c)));
        auto const s2 = 2.0 * x(c) * (y(a) - y(c)) / std::pow(x(a) - x(c), 2)
                        - (x(a) + x(c)) / (x(a) - x(c)) / x(c);
        auto const s3 = 2.0 * (y(a) - y(c)) / std::pow(x(a) - x(c), 3)
                        - (1.0 / x(a) + 1.0 / x(c)) / std::pow(x(a) - x(c), 2);
        auto const s4 = x(c) * s3;
        auto const s5 = x(c) * s3;
        auto const s6 = std::pow(x(c), 2) * s3;

        return s1 * dX2_dX(Be_trial) - s2 * Isym - s3 * outer_product(Be_trial)
               + s4 * outer_product(Be_trial, matrix3::Identity())
               + s5 * outer_product(matrix3::Identity(), Be_trial)
               + s6 * outer_product(matrix3::Identity());
    }
    // Derivative with all repeated eigenvalues
    std::cout << "Returning the symmetrix identity tensor and the first eigenvalue " << x(0) << "\n";
    return x.norm() < 1.0e-5 ? Isym : Isym / x(0);
}

matrix6 finite_strain_J2_plasticity::compute_B(matrix3 const& Be_trial) const
{
    matrix6 B(6, 6);
    // clang-format off
    B << 2 * Be_trial(0, 0),                  0,                  0,                  0, 2 * Be_trial(0, 2), 2 * Be_trial(0, 1),
                          0, 2 * Be_trial(1, 1),                  0, 2 * Be_trial(1, 2),                  0,                  0,
                          0,                  0, 2 * Be_trial(2, 2),                  0,                  0,                  0,
                          0, 2 * Be_trial(1, 2),                  0,     Be_trial(2, 2),                  0,                  0,
         2 * Be_trial(0, 2),                  0,                  0,                  0,     Be_trial(2, 2),      Be_trial(2, 1),
         2 * Be_trial(0, 1),                  0,                  0,                  0,     Be_trial(2, 1),      Be_trial(1, 1);
    // clang-format on
    return B;
}

matrix6 finite_strain_J2_plasticity::compute_H(matrix3 const& s) const
{
    matrix6 H(6, 6);
    // clang-format off
    H << 2 * s(0, 0),           0,           0,           0, 2 * s(0, 2), 2 * s(0, 1),
                   0, 2 * s(1, 1),           0, 2 * s(1, 2),           0,           0,
                   0,           0, 2 * s(2, 2),           0,           0,           0,
                   0, 2 * s(1, 2),           0,     s(2, 2),           0,           0,
         2 * s(0, 2),           0,           0,           0,     s(2, 2),     s(2, 1),
         2 * s(0, 1),           0,           0,           0,     s(1, 2),     s(2, 2);
    // clang-format on
    return H;
}

matrix6 finite_strain_J2_plasticity::dX2_dX(matrix3 const& X) const
{
    matrix6 dx2_dx(6, 6);
    // clang-format off
    dx2_dx << 2 * X(0, 0),         0.0,         0.0,               0.0,           X(2, 0),          X(1, 0),
                      0.0, 2 * X(1, 1),         0.0,           X(2, 1),               0.0,          X(1, 0),
                      0.0,         0.0, 2 * X(2, 2),           X(2, 1),           X(2, 0),              0.0,
                      0.0,     X(2, 1),     X(2, 1), X(1, 1) + X(2, 2),           X(1, 0),              0.0,
                  X(2, 0),           0,     X(2, 0),           X(1, 0), X(0, 0) + X(2, 2),           X(1, 2),
                  X(1, 0),     X(1, 0),         0.0,               0.0,           X(1, 2), X(0, 0) + X(1, 1);
    // clang-format on
    return dx2_dx;
}

matrix6 finite_strain_J2_plasticity::mandel_transformation() const
{
    matrix6 M = matrix6::Ones();

    M.block<3, 3>(0, 3) *= std::sqrt(2);
    M.block<3, 3>(3, 0) *= std::sqrt(2);
    M.block<3, 3>(3, 3) *= 2.0;

    return M;
}
}
