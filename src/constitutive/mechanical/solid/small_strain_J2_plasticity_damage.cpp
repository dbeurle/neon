#include "small_strain_J2_plasticity_damage.hpp"

#include "constitutive/internal_variables.hpp"
#include "numeric/mechanics"

#include "io/json.hpp"
#include "Exceptions.hpp"

#include <range/v3/view/transform.hpp>
#include <tbb/tbb.h>

#include <iostream>

namespace neon::mechanical::solid
{
small_strain_J2_plasticity_damage::small_strain_J2_plasticity_damage(
    std::shared_ptr<internal_variables_t>& variables, json const& material_data)
    : small_strain_J2_plasticity(variables, material_data), material(material_data)
{
    variables->add(internal_variables_t::Tensor::BackStress,
                   internal_variables_t::Tensor::KinematicHardening,
                   internal_variables_t::Scalar::Damage,
                   internal_variables_t::Scalar::EnergyReleaseRate);
}

small_strain_J2_plasticity_damage::~small_strain_J2_plasticity_damage() = default;

void small_strain_J2_plasticity_damage::update_internal_variables(double const time_step_size)
{
    // Retrieve the internal variables
    auto [plastic_strains,
          strains,
          cauchy_stresses,
          back_stresses,
          accumulated_kinematic_stresses] = variables
                                                ->fetch(internal_variables_t::Tensor::LinearisedPlasticStrain,
                                                        internal_variables_t::Tensor::LinearisedStrain,
                                                        internal_variables_t::Tensor::Cauchy,
                                                        internal_variables_t::Tensor::BackStress,
                                                        internal_variables_t::Tensor::KinematicHardening);

    // Retrieve the accumulated internal variables
    auto [accumulated_plastic_strains,
          von_mises_stresses,
          scalar_damages,
          energy_release_rates] = variables->fetch(internal_variables_t::Scalar::EffectivePlasticStrain,
                                                   internal_variables_t::Scalar::VonMisesStress,
                                                   internal_variables_t::Scalar::Damage,
                                                   internal_variables_t::Scalar::EnergyReleaseRate);

    auto& tangent_operators = variables->fetch(internal_variables_t::rank4::tangent_operator);

    // Compute the linear strain gradient from the displacement gradient
    strains = variables->fetch(internal_variables_t::Tensor::DisplacementGradient)
              | ranges::view::transform([](auto const& H) { return 0.5 * (H + H.transpose()); });

    // Perform the update algorithm for each quadrature point
    tbb::parallel_for(std::size_t{0}, strains.size(), [&](auto const l) {
        auto const& strain = strains[l];
        auto& plastic_strain = plastic_strains[l];
        auto& cauchy_stress = cauchy_stresses[l];
        auto& von_mises = von_mises_stresses[l];
        auto& back_stress = back_stresses[l];
        auto& scalar_damage = scalar_damages[l];
        auto& energy_var = energy_release_rates[l];

        // Elastic stress predictor
        cauchy_stress = compute_cauchy_stress(material.shear_modulus(),
                                              material.lambda(),
                                              strain - plastic_strain);

        matrix3 tau = deviatoric(cauchy_stress) / (1.0 - scalar_damage) - deviatoric(back_stress);

        // Trial von Mises stress
        von_mises = von_mises_stress(tau);

        energy_var = 0.5
                     * double_dot(strain - plastic_strain,
                                  compute_stress_like_matrix((1.0 - scalar_damage) * C_e,
                                                             strain - plastic_strain));

        // If this quadrature point is elastic, then set the tangent to the
        // elastic modulus and continue to the next quadrature point
        if (evaluate_yield_function(von_mises, back_stress) <= 0.0
            && evaluate_damage_yield_function(energy_var) <= 0.0)
        {
            tangent_operators[l] = C_e;
            return;
        }

        auto const plastic_increment = perform_radial_return(cauchy_stress,
                                                             back_stress,
                                                             scalar_damage,
                                                             accumulated_kinematic_stresses[l],
                                                             energy_var,
                                                             tangent_operators[l],
                                                             time_step_size,
                                                             strain - plastic_strain);

        tau = deviatoric(cauchy_stress) / (1.0 - scalar_damage) - deviatoric(back_stress);

        von_mises = von_mises_stress(tau);

        plastic_strain += plastic_increment * 3.0 / 2.0 * tau / (von_mises * (1 - scalar_damage));

        accumulated_plastic_strains[l] += plastic_increment / (1 - scalar_damage);

        // cauchy_stress, back_stress, scalar_damage, kin_hard, energy_var, C_algorithmic, are
        // updated within the radial_return routine std::cout << "delta_t  " << time_step_size <<
        // "\n";
    });
}

double small_strain_J2_plasticity_damage::perform_radial_return(matrix3& cauchy_stress,
                                                                matrix3& back_stress,
                                                                double& scalar_damage,
                                                                matrix3& kin_hard,
                                                                double& energy_var,
                                                                matrix6& tangent_operator,
                                                                double const& delta_t,
                                                                matrix3 const& eps_e_t)
{
    // TODO: initial guess could be computed base on a frozen yield surface (at
    // the increment n) to obtain good convergence. check computational methods for plasticity book

    vector16 y;
    y << 0.0, 0.0, voigt::kinetic::to(cauchy_stress), voigt::kinetic::to(back_stress),
        scalar_damage, energy_var;

    // The residual
    vector16 f = vector16::Ones();

    int iterations = 0, max_iterations = 25;

    matrix16 M = matrix16::Zero();

    M(0, 0) = M(1, 1) = M(14, 14) = M(15, 15) = 1.0;

    // TODO: double check kinetic kinematic
    M.block<6, 6>(2, 2) = voigt::kinetic::fourth_order_identity();

    auto const kp = material.plasticity_viscous_multiplier();
    auto const np = material.plasticity_viscous_exponent();
    auto const kd = material.damage_viscous_multiplier();
    auto const nd = material.damage_viscous_exponent();
    auto const C = material.kinematic_hardening_modulus();
    auto const gamma = material.softening_multiplier();

    // TODO: use relative error instead of absolute one (for f and delta_y)
    // TODO: use ResidualTolerance
    while (f.norm() > 1.0e-6 && iterations < max_iterations)
    {
        auto const d_lam_plastic = y(0);
        auto const d_lam_damage = y(1);

        auto const sigma = voigt::kinetic::from(y.segment<6>(2));
        auto const beta = voigt::kinetic::from(y.segment<6>(8));

        auto const d = y(14);
        auto const Y = y(15);

        matrix3 const tau = deviatoric(sigma) / (1.0 - d) - deviatoric(beta);

        auto const von_mises = von_mises_stress(tau);

        auto const f_vp = std::max(0.0, evaluate_yield_function(von_mises, beta));
        auto const f_d = std::max(0.0, evaluate_damage_yield_function(Y));

        matrix3 const normal_tild = 3.0 / 2.0 * tau / von_mises;
        matrix3 const normal = normal_tild / (1.0 - d);

        // The residual TODO: compute_stress_like_vector and compute_stress_like_matrix could be
        // optimised
        // clang-format off
        f << d_lam_plastic - kp * std::pow(f_vp, np),
             d_lam_damage - kd * std::pow(f_d, nd),
             voigt::kinetic::to(sigma) - (1.0 - d) * compute_stress_like_vector(C_e, eps_e_t) + delta_t * d_lam_plastic * compute_stress_like_vector(C_e, normal_tild),
             voigt::kinetic::to(beta - C * kin_hard - C * delta_t * d_lam_plastic * (normal_tild - gamma / C * beta)),
             d - scalar_damage - delta_t * d_lam_damage,
             Y - 0.5   * double_dot(eps_e_t - delta_t * d_lam_plastic * normal_tild / (1.0 - d),
             compute_stress_like_matrix(C_e, eps_e_t - delta_t * d_lam_plastic * normal_tild / (1.0 - d)));

        M.block<1, 6>(0, 2) = voigt::kinetic::to(deviatoric(-np * kp * std::pow(f_vp, np - 1) * normal));
        M.block<1, 6>(0, 8) = voigt::kinetic::to(np * kp * std::pow(f_vp, np - 1) * (normal_tild - gamma / C * beta));

        M(0, 14) = -np * kp * std::pow(f_vp, np - 1) * double_dot(normal_tild, deviatoric(sigma)) / std::pow(1.0 - d, 2);
        M(1, 15) = -nd * kd * std::pow(f_d, nd - 1);

        M.block<6, 1>(2, 0) = delta_t * compute_stress_like_vector(C_e, normal_tild);
        M.block<6, 1>(2, 14) = compute_stress_like_vector(C_e, eps_e_t);
        M.block<6, 1>(8, 0) = voigt::kinetic::to(-delta_t * C * (normal_tild - gamma / C * beta));
        M.block<6, 6>(8, 8) = voigt::kinetic::fourth_order_identity() * (1.0 + delta_t * d_lam_plastic * gamma);

        M(14, 1) = -delta_t;

        M(15, 0) = double_dot(delta_t * normal_tild / (1.0 - d),
                              compute_stress_like_matrix(C_e, eps_e_t - delta_t * d_lam_plastic * normal_tild / (1.0 - d)));

        M(15, 14) = delta_t * d_lam_plastic / std::pow(1.0 - d, 2)
                    * double_dot(normal_tild,
                                 compute_stress_like_matrix(C_e, eps_e_t - delta_t * d_lam_plastic * normal_tild / (1.0 - d)));
        // clang-format on

        // Eigen::JacobiSVD<Eigen::MatrixXd> svd(M);
        // double cond = svd.singularValues()(0) /
        // svd.singularValues()(svd.singularValues().size()
        // - 1); std::cout << cond << "\n\n";
        // vector a = -M.inverse() * f;

        y -= M.fullPivLu().solve(f);

        iterations++;
    }

    if (iterations == max_iterations)
    {
        std::cout << "MAXIMUM NUMBER OF ITERATIONS IN RADIAL RETURN REACHED\n";
        std::cout << "The residual norm " << f.norm() << "\n";
        // std::cout << "The increment norm " << a.norm() << "\n";
        throw computational_error("Radial return failure\n");
    }

    cauchy_stress = voigt::kinetic::from(y.segment<6>(2));

    back_stress = voigt::kinetic::from(y.segment<6>(8));

    scalar_damage = y(14);

    energy_var = y(15);

    kin_hard = back_stress / C;

    tangent_operator = (1.0 - scalar_damage) * mandel_notation(C_e)
                       * mandel_notation(M.inverse().block<6, 6>(2, 2));

    return y(0) * delta_t; // plastic_increment
}

double small_strain_J2_plasticity_damage::evaluate_yield_function(double const von_mises,
                                                                  matrix3 const& back_stress) const
{
    return von_mises
           + 0.5 * material.softening_multiplier() / material.kinematic_hardening_modulus()
                 * double_dot(back_stress, back_stress)
           - material.yield_stress(0.0);
}

double small_strain_J2_plasticity_damage::evaluate_damage_yield_function(double const energy_var) const
{
    return energy_var - std::pow(material.yield_stress(0.0), 2) / (2 * material.elastic_modulus());
}

matrix3 small_strain_J2_plasticity_damage::compute_stress_like_matrix(matrix6 const& tangent_operator,
                                                                      matrix3 const& strain_like) const
{
    // could get the tangent_operator directly without passing it to this function
    return voigt::kinetic::from(tangent_operator * voigt::kinematic::to(strain_like));
}

vector6 small_strain_J2_plasticity_damage::compute_stress_like_vector(matrix6 const& tangent_operator,
                                                                      matrix3 const& strain_like) const
{
    // could get the tangent_operator directly without passing it to this function
    return tangent_operator * voigt::kinematic::to(strain_like);
}
}
