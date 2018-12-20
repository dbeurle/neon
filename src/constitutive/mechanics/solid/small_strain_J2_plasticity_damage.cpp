#include "small_strain_J2_plasticity_damage.hpp"

#include "constitutive/internal_variables.hpp"
#include "numeric/mechanics"

#include "io/json.hpp"
#include "exceptions.hpp"

#include <range/v3/view/transform.hpp>
#include <tbb/parallel_for.h>

#include <iostream>

namespace neon::mechanics::solid
{
small_strain_J2_plasticity_damage::small_strain_J2_plasticity_damage(
    std::shared_ptr<internal_variables_t>& variables,
    json const& material_data)
    : small_strain_J2_plasticity(variables, material_data), material(material_data)
{
    variables->add(variable::second::back_stress,
                   variable::second::kinematic_hardening,
                   variable::scalar::damage,
                   variable::scalar::energy_release_rate);

    names.emplace("back_stress");
    names.emplace("kinematic_hardening");
    names.emplace("damage");
    names.emplace("energy_release_rate");
}

small_strain_J2_plasticity_damage::~small_strain_J2_plasticity_damage() = default;

void small_strain_J2_plasticity_damage::update_internal_variables(double const time_step_size)
{
    // Retrieve the internal variables
    auto& plastic_strains = variables->get(variable::second::linearised_plastic_strain);
    auto& strains = variables->get(variable::second::linearised_strain);
    auto& accumulated_plastic_strains = variables->get(variable::scalar::effective_plastic_strain);
    auto& kinematic_hardening = variables->get(variable::second::kinematic_hardening);

    auto& cauchy_stresses = variables->get(variable::second::cauchy_stress);
    auto& back_stresses = variables->get(variable::second::back_stress);
    auto& von_mises_stresses = variables->get(variable::scalar::von_mises_stress);

    auto& energy_release_rates = variables->get(variable::scalar::energy_release_rate);
    auto& scalar_damages = variables->get(variable::scalar::damage);

    auto& tangent_operators = variables->get(variable::fourth::tangent_operator);

    // Compute the linear strain gradient from the displacement gradient
    strains = variables->get(variable::second::displacement_gradient)
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
        cauchy_stress = compute_cauchy_stress_scalar_damage(material.shear_modulus(),
                                                            material.lambda(),
                                                            strain - plastic_strain,
                                                            scalar_damage);

        // Trial state
        matrix3 relative_stress = cauchy_stress - back_stress;

        von_mises = von_mises_stress(relative_stress);

        energy_var = compute_energy_release_rate(C_e, strain - plastic_strain);

        // If this quadrature point is elastic, then set the tangent to the
        // elastic modulus and continue to the next quadrature point
        if (evaluate_yield_function(von_mises, back_stress, scalar_damage) <= 0.0
            && evaluate_damage_yield_function(energy_var) <= 0.0)
        {
            tangent_operators[l] = C_e;
            return;
        }

        // TODO: here the normal is considered to be non-constant. However, only its magnitude is
        // changing so only its derivative w.r.t damage is non-zero
        auto const plastic_increment = perform_radial_return(cauchy_stress,
                                                             back_stress,
                                                             scalar_damage,
                                                             kinematic_hardening[l],
                                                             energy_var,
                                                             tangent_operators[l],
                                                             time_step_size,
                                                             strain - plastic_strain);

        relative_stress = cauchy_stress - back_stress;

        von_mises = von_mises_stress(relative_stress);

        // Compute the normal direction to the yield surface
        matrix3 const normal = std::sqrt(3.0 / 2.0) * deviatoric(relative_stress)
                               / (deviatoric(relative_stress).norm() * (1 - scalar_damage));

        plastic_strain += plastic_increment * normal;

        accumulated_plastic_strains[l] += plastic_increment / (1 - scalar_damage);
    });
}

double small_strain_J2_plasticity_damage::perform_radial_return(matrix3& cauchy_stress_t,
                                                                matrix3& back_stress_t,
                                                                double& scalar_damage_t,
                                                                matrix3& kin_hard_t,
                                                                double& energy_var_t,
                                                                matrix6& tangent_operator,
                                                                double const& delta_t,
                                                                matrix3 const& eps_e_t)
{
    // TODO: initial guess could be computed base on a frozen yield surface (at
    // the preceding time increment) to obtain good convergence. check "computational methods for plasticity" book

    vector16 solution;
    solution << 0.0, 0.0, voigt::kinetic::to(cauchy_stress_t), voigt::kinetic::to(back_stress_t),
        scalar_damage_t, energy_var_t;

    vector16 residual = vector16::Ones();

    matrix16 jacobian = matrix16::Zero();
    jacobian(0, 0) = jacobian(1, 1) = jacobian(14, 14) = jacobian(15, 15) = 1.0;
    jacobian.block<6, 6>(2, 2) = voigt::kinetic::fourth_order_identity();
    jacobian.block<6, 6>(8, 8) = voigt::kinetic::fourth_order_identity();

    auto const s_p = material.plasticity_viscous_denominator();
    auto const n_p = material.plasticity_viscous_exponent();
    auto const k_p = std::pow(s_p, -n_p); // plasticity viscous multiplier
    auto const s_d = material.damage_viscous_denominator();
    auto const n_d = material.damage_viscous_exponent();
    auto const k_d = std::pow(s_d, -n_d); // damage viscous multiplier
    auto const c_p = material.kinematic_hardening_modulus();
    auto const a_p = material.softening_multiplier();

    // TODO: use relative error instead of absolute one (for residual and delta_y)
    // Perform the non-linear hardening solve
    int iterations = 0;
    auto constexpr max_iterations{50};
    while (residual.norm() > 1.0e-6 && iterations < max_iterations)
    {
        auto const delta_plastic_multiplier = solution(0);
        auto const delta_damage_multiplier = solution(1);

        auto const cauchy_stress = voigt::kinetic::from(solution.segment<6>(2));
        auto const back_stress = voigt::kinetic::from(solution.segment<6>(8));

        auto const scalar_damage = solution(14);
        auto const one_minus_damage = (1.0 - scalar_damage);
        auto const energy = solution(15);

        // Compute the normal direction to the yield surface
        matrix3 const relative_stress = cauchy_stress - back_stress;
        auto const von_mises = von_mises_stress(relative_stress);
        matrix3 const normal = std::sqrt(3.0 / 2.0) * deviatoric(relative_stress)
                               / (deviatoric(relative_stress).norm() * one_minus_damage);

        matrix3 const eps_e = eps_e_t - delta_plastic_multiplier * normal;

        auto const yield_function = std::max(0.0,
                                             evaluate_yield_function(von_mises,
                                                                     back_stress,
                                                                     scalar_damage));
        auto const damage_yield_function = std::max(0.0, evaluate_damage_yield_function(energy));

        matrix3 const& d_yield_d_stress = normal;
        matrix3 const d_yield_d_back_stress = a_p / c_p * back_stress - normal;
        double const d_yield_d_damage = von_mises / std::pow(one_minus_damage, 2);

        // clang-format off
        matrix6 const d_normal_wrt_cauchy_stress = std::sqrt(3.0 / 2.0) / one_minus_damage
                                                 * (voigt::kinetic::deviatoric() / deviatoric(relative_stress).norm()
                                                    - outer_product(deviatoric(relative_stress), deviatoric(relative_stress))
                                                    / std::pow(deviatoric(relative_stress).norm(), 3));
        matrix6 const d_normal_wrt_back_stress = -d_normal_wrt_cauchy_stress;
        matrix3 const d_normal_wrt_damage = normal / one_minus_damage;

        residual << delta_plastic_multiplier - k_p * std::pow(yield_function, n_p) * delta_t,
                    delta_damage_multiplier - k_d * std::pow(damage_yield_function, n_d) * delta_t,
                    voigt::kinetic::to(cauchy_stress - compute_cauchy_stress_scalar_damage(
                                        material.shear_modulus(),material.lambda(), eps_e, scalar_damage)),
                    voigt::kinetic::to(back_stress - back_stress_t
                                      - delta_plastic_multiplier * (c_p * normal - a_p * back_stress)),
                    scalar_damage - scalar_damage_t - delta_damage_multiplier,
                    energy - compute_energy_release_rate(C_e, eps_e);

        // Derivative of the plastic multiplier equation
        jacobian.block<1, 6>(0, 2) = voigt::kinetic::to(
            -n_p * k_p * std::pow(yield_function, n_p - 1) * d_yield_d_stress * delta_t);
        jacobian.block<1, 6>(0, 8) = voigt::kinetic::to(
            -n_p * k_p * std::pow(yield_function, n_p - 1) * d_yield_d_back_stress * delta_t);
        jacobian(0, 14) =
            -n_p * k_p * std::pow(yield_function, n_p - 1) * d_yield_d_damage * delta_t;

        // Derivative of the damage multiplier equation
        jacobian(1, 15) = -n_d * k_d * std::pow(damage_yield_function, n_d - 1) * delta_t;

        // Derivative of the stress equation
        jacobian.block<6, 1>(2, 0) = compute_stress_like_vector_voigt(one_minus_damage * C_e, normal);
        jacobian.block<6, 6>(2, 2) = voigt::kinetic::fourth_order_identity()
                                     + mandel_notation(one_minus_damage * C_e)
                                     * mandel_notation(delta_plastic_multiplier * d_normal_wrt_cauchy_stress);
        jacobian.block<6, 6>(2, 8) = mandel_notation(one_minus_damage * C_e)
                                     * mandel_notation(delta_plastic_multiplier * d_normal_wrt_back_stress);
        jacobian.block<6, 1>(2, 14) = compute_stress_like_vector_voigt(C_e,
                                          one_minus_damage * delta_plastic_multiplier * d_normal_wrt_damage + eps_e);

        // Derivative of the backstress equation
        jacobian.block<6, 1>(8, 0) = voigt::kinetic::to(a_p * back_stress - c_p * normal);
        jacobian.block<6, 6>(8, 2) = -c_p * delta_plastic_multiplier * d_normal_wrt_cauchy_stress;
        jacobian.block<6, 6>(8, 8) = voigt::kinetic::fourth_order_identity()
                                     * (1.0 + a_p * delta_plastic_multiplier)
                                     - c_p * delta_plastic_multiplier * d_normal_wrt_back_stress;
        jacobian.block<6, 1>(8, 14) = voigt::kinetic::to(-c_p * delta_plastic_multiplier * d_normal_wrt_damage);

        // Derivative of the damage equation
        jacobian(14, 1) = -1;

        // Derivative of the energy release rate equation
        jacobian(15, 0) = double_dot(normal,
                                 compute_stress_like_matrix(C_e, eps_e));
        jacobian.block<1, 6>(15, 2) = delta_plastic_multiplier
                                      * (d_normal_wrt_cauchy_stress
                                      * compute_stress_like_vector_voigt(C_e, eps_e)).transpose();
        jacobian.block<1, 6>(15, 8) = delta_plastic_multiplier
                                      * (d_normal_wrt_back_stress
                                      * compute_stress_like_vector_voigt(C_e, eps_e)).transpose();
        jacobian(15, 14) = delta_plastic_multiplier
                           * double_dot(d_normal_wrt_damage, compute_stress_like_matrix(C_e, eps_e));

        // clang-format on

        solution -= jacobian.fullPivLu().solve(residual);

        iterations++;
    }

    // ensure that the damage does not exceed one
    if (solution(14) > 1)
    {
        std::cout << "Damage is greater than one\nConsider decreasing the load\n";
        throw computational_error("Radial return failure\n");
    }

    if (iterations == max_iterations)
    {
        std::cout << "MAXIMUM NUMBER OF ITERATIONS IN RADIAL RETURN IS REACHED\n";
        std::cout << "The residual norm " << residual.norm() << "\n";
        throw computational_error("Radial return failure\n");
    }

    // ensure a positive value for the plastic multiplier
    auto plastic_increment = solution(0) > 0.0 ? solution(0) : 0.0;

    cauchy_stress_t = voigt::kinetic::from(solution.segment<6>(2));

    back_stress_t = voigt::kinetic::from(solution.segment<6>(8));

    // ensure no damage-reduction/healing
    scalar_damage_t = solution(14) > scalar_damage_t ? solution(14) : scalar_damage_t;

    energy_var_t = solution(15);

    kin_hard_t = back_stress_t / c_p;

    tangent_operator = (1.0 - scalar_damage_t) * mandel_notation(C_e)
                       * mandel_notation(jacobian.inverse().block<6, 6>(2, 2));

    return plastic_increment * delta_t; // plastic_increment
}

double small_strain_J2_plasticity_damage::evaluate_yield_function(double const von_mises,
                                                                  matrix3 const& back_stress,
                                                                  double const damage) const
{
    return von_mises / (1 - damage)
           + 0.5 * material.softening_multiplier() / material.kinematic_hardening_modulus()
                 * double_dot(back_stress, back_stress)
           - material.yield_stress(0.0);
}

double small_strain_J2_plasticity_damage::evaluate_damage_yield_function(double const energy_var) const
{
    return energy_var - std::pow(material.yield_stress(0.0), 2) / (2 * material.elastic_modulus());
}

double small_strain_J2_plasticity_damage::compute_energy_release_rate(matrix6 const& tangent_operator,
                                                                      matrix3 const& elastic_strain) const
{
    return 0.5
           * double_dot(elastic_strain, compute_stress_like_matrix(tangent_operator, elastic_strain));
}

matrix3 small_strain_J2_plasticity_damage::compute_stress_like_matrix(matrix6 const& tangent_operator,
                                                                      matrix3 const& strain_like) const
{
    return voigt::kinetic::from(compute_stress_like_vector_voigt(tangent_operator, strain_like));
}

vector6 small_strain_J2_plasticity_damage::compute_stress_like_vector_voigt(
    matrix6 const& tangent_operator,
    matrix3 const& strain_like) const
{
    return tangent_operator * voigt::kinematic::to(strain_like);
}
}
