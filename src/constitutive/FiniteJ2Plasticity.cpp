
#include "FiniteJ2Plasticity.hpp"

#include "Exceptions.hpp"
#include "InternalVariables.hpp"

#include <iostream>
#include <json/value.h>
#include <range/v3/view.hpp>
#include <unsupported/Eigen/MatrixFunctions>

namespace neon
{
FiniteJ2Plasticity::FiniteJ2Plasticity(InternalVariables& variables, Json::Value const& material_data)
    : J2Plasticity(variables, material_data)
{
    variables.add(InternalVariables::Scalar::VonMisesStress,
                  InternalVariables::Scalar::EffectivePlasticStrain);

    variables.add(InternalVariables::Tensor::HenckyStrainElastic);

    // Add material tangent with the linear elasticity moduli
    variables.add(InternalVariables::Matrix::TruesdellModuli,
                  consistent_tangent(1.0, Matrix3::Zero(), Matrix3::Zero(), C_e));

    // std::cout << "Constructed finite strain J2 plasticity model\n";
}

FiniteJ2Plasticity::~FiniteJ2Plasticity() = default;

void FiniteJ2Plasticity::update_internal_variables(double const time_step_size)
{
    using namespace ranges;

    auto const shear_modulus = material.shear_modulus();

    // Extract the internal variables
    auto[F_list,
         log_strain_e_list,
         cauchy_stress_list] = variables(InternalVariables::Tensor::DeformationGradient,
                                         InternalVariables::Tensor::HenckyStrainElastic,
                                         InternalVariables::Tensor::Cauchy);

    auto const F_old_list = variables[InternalVariables::Tensor::DeformationGradient];

    auto const J_list = variables(InternalVariables::Scalar::DetF);

    // Retrieve the accumulated internal variables
    auto[accumulated_plastic_strain_list,
         von_mises_list] = variables(InternalVariables::Scalar::EffectivePlasticStrain,
                                     InternalVariables::Scalar::VonMisesStress);

    auto& C_list = variables(InternalVariables::Matrix::TruesdellModuli);

    auto const F_inc_list = view::zip(F_list, F_old_list) | view::transform([](auto const& tpl) {
                                auto const & [ F, F_old ] = tpl;
                                return F * F_old.inverse();
                            });

    // Perform the update algorithm for each quadrature point
    for (auto l = 0; l < F_list.size(); l++)
    {
        auto const& F_inc = F_inc_list[l];
        auto const J = J_list[l];

        auto& cauchy_stress = cauchy_stress_list[l];
        auto& accumulated_plastic_strain = accumulated_plastic_strain_list[l];
        auto& von_mises = von_mises_list[l];
        auto& log_strain_e = log_strain_e_list[l];

        // std::cout << "log strain elastic\n" << log_strain_e << std::endl;

        // Elastic trial deformation gradient
        Matrix3 const B_e = (2.0 * log_strain_e).exp();

        // std::cout << "B elastic\n" << B_e << std::endl;

        // std::cout << "Finc\n" << F_inc << std::endl;

        // Elastic trial left Cauchy-Green deformation tensor
        Matrix3 const B_e_trial = F_inc * B_e * F_inc.transpose();

        // std::cout << "B elastic trial\n" << B_e_trial << std::endl;

        // Trial Logarithmic elastic strain
        log_strain_e = 0.5 * B_e_trial.log();

        // std::cout << "log strain elastic\n" << log_strain_e << std::endl;

        // Elastic stress predictor
        cauchy_stress = compute_cauchy_stress(log_strain_e) / J;

        // Trial von Mises stress
        von_mises = von_mises_stress(cauchy_stress);

        // std::cout << "von Mises " << von_mises << std::endl;

        // Compute the initial estimate of the yield function for the material
        // and decide if the stress return needs to be computed
        if (auto const f = evaluate_yield_function(von_mises, accumulated_plastic_strain); f <= 0.0)
        {
            // std::cout << "f = " << f << std::endl;
            C_list[l] = consistent_tangent(J, log_strain_e, cauchy_stress, C_e);
            continue;
        }

        auto const von_mises_trial = von_mises;

        // std::cout << "\nQUADRATURE POINT PLASTIC\n";

        // Compute the normal direction to the yield surface which remains
        // constant throughout the radial return method
        Matrix3 const normal = deviatoric(cauchy_stress) / deviatoric(cauchy_stress).norm();

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
        auto const D_ep = algorithmic_tangent(plastic_increment,
                                              accumulated_plastic_strain,
                                              von_mises_trial,
                                              normal);

        // Compute the elastic-plastic tangent modulus for large strain
        C_list[l] = consistent_tangent(J, log_strain_e, cauchy_stress, D_ep);
    }
}

CMatrix FiniteJ2Plasticity::consistent_tangent(double const J,
                                               Matrix3 const& Be_trial,
                                               Matrix3 const& cauchy_stress,
                                               CMatrix const& C)
{
    // Convert to Mandel notation so matrix multiplication == double dot operation
    CMatrix const D = voigt_to_mandel(C);
    CMatrix const L = voigt_to_mandel(compute_L(Be_trial));
    CMatrix const B = voigt_to_mandel(compute_B(Be_trial));

    CMatrix const H = compute_H(cauchy_stress);

    return 1.0 / (2.0 * J) * D * L * B - H;
}

std::tuple<Vector3, std::array<Matrix3, 3>, bool, std::array<int, 3>> FiniteJ2Plasticity::
    compute_eigenvalues_eigenprojections(Matrix3 const& X) const
{
    Eigen::SelfAdjointEigenSolver<Matrix3> eigen_solver(X);

    if (eigen_solver.info() != Eigen::Success)
    {
        throw computational_error("Eigenvalue solver failed in finite plasticity routine\n");
    }

    Vector3 const x = eigen_solver.eigenvalues();
    Matrix3 const v = eigen_solver.eigenvectors();

    // Eigenprojections
    std::array<Matrix3, 3> E = {{v.col(0) * v.col(0).transpose(),
                                 v.col(1) * v.col(1).transpose(),
                                 v.col(2) * v.col(2).transpose()}};

    bool is_repeated = false;
    std::array<int, 3> abc_ordering{{0, 1, 2}};

    // Need to check if there are repeated eigenvalues to machine precision
    // if (x.norm() < 1.0e-2 || (is_approx(x(0), x(1)) && is_approx(x(1), x(2))))
    if (is_approx(x(0), x(1)) && is_approx(x(1), x(2)))
    {
        std::cout << "We have repeated eigenvalues!\n" << x << std::endl;
        E[0] = Matrix3::Identity();
        is_repeated = true;
        abc_ordering = {{-1, -1, -1}};
    }
    else if (!is_approx(x(0), x(1)) && is_approx(x(1), x(2)))
    {
        E[1] = Matrix3::Identity() - E[0];
        is_repeated = true;
        abc_ordering = {{0, 1, 2}};
    }
    else if (!is_approx(x(2), x(0)) && is_approx(x(0), x(1)))
    {
        E[0] = Matrix3::Identity() - E[2];
        is_repeated = true;
        abc_ordering = {{2, 0, 1}};
    }
    else if (!is_approx(x(1), x(2)) && is_approx(x(2), x(0)))
    {
        E[2] = Matrix3::Identity() - E[1];
        is_repeated = true;
        abc_ordering = {{1, 2, 0}};
    }
    return {x, E, is_repeated, abc_ordering};
}

CMatrix FiniteJ2Plasticity::derivative_tensor_log_unique(Matrix3 const& Be_trial,
                                                         Vector3 const& x,
                                                         std::array<Matrix3, 3> const& E,
                                                         std::array<int, 3> const& abc_ordering) const
{
    auto const[a, b, c] = abc_ordering;
    return std::log(x(a)) / ((x(a) - x(b)) * (x(a) - x(c)))
               * (dX2_dX(Be_trial) - (x(b) - x(c)) * Isym
                  - (2.0 * x(a) - x(b) - x(c)) * outer_product(E[a])
                  - (x(b) - x(c)) * (outer_product(E[b]) - outer_product(E[c])))
           + 1.0 / x(a) * outer_product(E[a]);
}

CMatrix FiniteJ2Plasticity::compute_L(Matrix3 const& Be_trial) const
{
    auto const[x, E, is_repeated, abc_ordering] = compute_eigenvalues_eigenprojections(Be_trial);

    if (!is_repeated)
    {
        return derivative_tensor_log_unique(Be_trial, x, E, {{0, 1, 2}})
               + derivative_tensor_log_unique(Be_trial, x, E, {{2, 0, 1}})
               + derivative_tensor_log_unique(Be_trial, x, E, {{1, 2, 0}});
    }
    else if (is_repeated && abc_ordering[0] >= 0)
    {
        // Derivative when there is one repeated eigenvalue
        auto const & [ a, b, c ] = abc_ordering;

        std::cout << "Eigenvalues\n" << x << std::endl;

        Vector3 const y = x.array().log();

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
               + s4 * outer_product(Be_trial, Matrix3::Identity())
               + s5 * outer_product(Matrix3::Identity(), Be_trial)
               + s6 * outer_product(Matrix3::Identity());
    }
    // Derivative with all repeated eigenvalues
    std::cout << "Returning the symmetrix identity tensor and the first eigenvalue " << x(0) << "\n";
    return x.norm() < 1.0e-5 ? Isym : Isym / x(0);
}

CMatrix FiniteJ2Plasticity::compute_B(Matrix3 const& Be_trial) const
{
    CMatrix B(6, 6);
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

CMatrix FiniteJ2Plasticity::compute_H(Matrix3 const& s) const
{
    CMatrix H(6, 6);
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

CMatrix FiniteJ2Plasticity::dX2_dX(Matrix3 const& X) const
{
    CMatrix dx2_dx(6, 6);
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

CMatrix FiniteJ2Plasticity::mandel_transformation() const
{
    CMatrix M = CMatrix::Ones(6, 6);

    M.block<3, 3>(0, 3) *= std::sqrt(2);
    M.block<3, 3>(3, 0) *= std::sqrt(2);
    M.block<3, 3>(3, 3) *= 2.0;

    return M;
}
}

// auto const I1 = X.trace();
// auto const I2 = 0.5 * (std::pow(X.trace(), 2) - (X * X).trace());
// auto const I3 = X.determinant();
//
// auto const R = (-2.0 * I1 + 9.0 * I1 * I2 - 27.0 * I3) / 54.0;
// auto const Q = (std::pow(I1, 2) - 3.0 * I2) / 9.0;
//
// auto const theta = std::acos(R / std::sqrt(std::pow(Q, 3)));
//
// Vector3 x(-2.0 * std::sqrt(Q) * std::cos(theta / 3.0) + I1 / 3.0,
//           -2.0 * std::sqrt(Q) * std::cos((theta + 2.0 * M_PI) / 3.0) + I1 / 3.0,
//           -2.0 * std::sqrt(Q) * std::cos((theta - 2.0 * M_PI) / 3.0) + I1 / 3.0);
//
// std::array<Matrix3, 3> E;
// int repeated_eigenvalues = 0;
//
// // Perform checks to determine if there are repeated eigenvalues
// if (!is_approx(x1, x2) && !is_approx(x2, x3) && !is_approx(x1, x3))
// {
//     for (int i = 0; i < 3; i++)
//     {
//         E[i] = x(i) / (2.0 * std::pow(x(i), 3) - I1 * std::pow(x(i), 2) + I3)
//                * (X * X - (I1 - x(i)) * X + I3 / x(i) * Matrix3::Identity());
//     }
// }
// else if (!is_approx(x1, x2) && is_approx(x2, x3))
// {
//     repeated_eigenvalue = 1;
// }
// else
// {
//     E[0] = Matrix3::Identity();
//     repeated_eigenvalue = 2;
// }
