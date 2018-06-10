
#include "constitutive/mechanical/solid/gaussian_ageing_affine_microsphere.hpp"

#include "constitutive/internal_variables.hpp"
#include "constitutive/mechanical/detail/microsphere.hpp"
#include "constitutive/mechanical/volumetric_free_energy.hpp"
#include "quadrature/trapezoidal.hpp"
#include "solver/time/runge_kutta_integration.hpp"
#include "numeric/float_compare.hpp"

#include <tbb/parallel_for.h>

#include <cassert>
#include <numeric>

namespace neon::mechanical::solid
{
gaussian_ageing_affine_microsphere::gaussian_ageing_affine_microsphere(
    std::shared_ptr<internal_variables_t>& variables,
    json const& material_data,
    unit_sphere_quadrature::point const rule)
    : gaussian_affine_microsphere{variables, material_data, rule}, material{material_data}
{
    variables->add(internal_variables_t::scalar::active_shear_modulus,
                   internal_variables_t::scalar::inactive_shear_modulus,
                   internal_variables_t::scalar::active_segments,
                   internal_variables_t::scalar::reduction_factor,
                   internal_variables_t::scalar::inactive_segments);

    variables->add(internal_variables_t::vector::first_ageing_moduli,
                   internal_variables_t::vector::second_ageing_moduli,
                   internal_variables_t::vector::third_ageing_moduli,
                   internal_variables_t::vector::first_previous,
                   internal_variables_t::vector::second_previous,
                   internal_variables_t::vector::third_previous);

    for (auto& values : variables->get(internal_variables_t::vector::first_ageing_moduli))
    {
        values.resize(unit_sphere.points(), 0.0);
    }
    for (auto& values : variables->get(internal_variables_t::vector::second_ageing_moduli))
    {
        values.resize(unit_sphere.points(), 0.0);
    }
    for (auto& values : variables->get(internal_variables_t::vector::third_ageing_moduli))
    {
        values.resize(unit_sphere.points(), 0.0);
    }
    for (auto& values : variables->get(internal_variables_t::vector::first_previous))
    {
        values.resize(unit_sphere.points(), 0.0);
    }
    for (auto& values : variables->get(internal_variables_t::vector::second_previous))
    {
        values.resize(unit_sphere.points(), 0.0);
    }
    for (auto& values : variables->get(internal_variables_t::vector::third_previous))
    {
        values.resize(unit_sphere.points(), 0.0);
    }

    auto [active_shear_modulus,
          active_segments,
          reduction] = variables->get(internal_variables_t::scalar::active_shear_modulus,
                                      internal_variables_t::scalar::active_segments,
                                      internal_variables_t::scalar::reduction_factor);

    std::fill(begin(active_shear_modulus), end(active_shear_modulus), material.shear_modulus());
    std::fill(begin(active_segments), end(active_segments), material.segments_per_chain());
    std::fill(begin(reduction), end(reduction), 1.0);

    variables->commit();
}

void gaussian_ageing_affine_microsphere::update_internal_variables(double const time_step_size)
{
    auto [active_shear_modulus,
          inactive_shear_modulus,
          reduction_factor,
          active_segments,
          inactive_segments] = variables->get(internal_variables_t::scalar::active_shear_modulus,
                                              internal_variables_t::scalar::inactive_shear_modulus,
                                              internal_variables_t::scalar::reduction_factor,
                                              internal_variables_t::scalar::active_segments,
                                              internal_variables_t::scalar::inactive_segments);

    // accumulated integral values from time integration
    auto [moduli1,
          moduli2,
          moduli3] = variables->get(internal_variable_t::vector::first_ageing_moduli,
                                    internal_variable_t::vector::second_ageing_moduli,
                                    internal_variable_t::vector::third_ageing_moduli);

    auto [sphere1, sphere2, sphere3] = variables->get(internal_variable_t::vector::first_previous,
                                                      internal_variable_t::vector::second_previous,
                                                      internal_variable_t::vector::third_previous);

    auto& moduli1_old = variables->get_old(internal_variable_t::vector::first_ageing_moduli);
    auto& moduli2_old = variables->get_old(internal_variable_t::vector::second_ageing_moduli);
    auto& moduli3_old = variables->get_old(internal_variable_t::vector::third_ageing_moduli);

    auto& sphere1_old = variables->get_old(internal_variable_t::vector::first_previous);
    auto& sphere2_old = variables->get_old(internal_variable_t::vector::second_previous);
    auto& sphere3_old = variables->get_old(internal_variable_t::vector::third_previous);

    auto& cauchy_stresses = variables->get(internal_variables_t::second::cauchy_stress);

    auto const& det_F = variables->get(internal_variables_t::scalar::DetF);

    auto const& deformation_gradients = variables->get(
        internal_variables_t::second::DeformationGradient);

    auto& tangent_operators = variables->get(internal_variables_t::fourth::tangent_operator);

    auto const K{material.bulk_modulus()};

    tbb::parallel_for(std::size_t{0}, deformation_gradients.size(), [&, this](auto const l) {
        matrix3 const F_bar = unimodular(deformation_gradients[l]);

        auto& modulus1 = moduli1[l];
        auto& modulus2 = moduli2[l];
        auto& modulus3 = moduli3[l];

        auto const& modulus1_old = moduli1_old[l];
        auto const& modulus2_old = moduli2_old[l];
        auto const& modulus3_old = moduli3_old[l];

        if (!is_approx(time_step_size, 0.0))
        {
            // Integrate the network evolution differential equations through
            // the micromechanical material by passing our network parameters
            vector5 parameters;
            parameters(0) = active_shear_modulus[l];
            parameters(1) = inactive_shear_modulus[l];
            parameters(2) = reduction_factor[l];
            parameters(3) = active_segments[l];
            parameters(4) = inactive_segments[l];

            parameters = material.integrate(parameters, time_step_size);

            // Update the history variables for plotting
            active_shear_modulus[l] = parameters(0);
            inactive_shear_modulus[l] = parameters(1);
            reduction_factor[l] = parameters(2);
            active_segments[l] = parameters(3);
            inactive_segments[l] = parameters(4);

            last_time_step_size = time_step_size;
        }

        // Partially integrate secondary modulus
        unit_sphere.for_each([&](auto const& coordinates, auto const index) {
            auto const& [r, _] = coordinates;

            auto const micro_stretch = compute_microstretch(deformed_tangent(F_bar, r));

            auto const creation_rate = material.creation_rate(active_shear_modulus[l],
                                                              inactive_shear_modulus[l],
                                                              active_segments[l],
                                                              inactive_segments[l]);

            auto const sphere1_eval = evaluate_integrand_first(creation_rate,
                                                               reduction_factor[l],
                                                               micro_stretch);

            auto const sphere2_eval = evaluate_integrand_second(creation_rate,
                                                                reduction_factor[l],
                                                                micro_stretch);

            auto const sphere3_eval = evaluate_integrand_third(creation_rate,
                                                               reduction_factor[l],
                                                               micro_stretch);

            modulus1.at(index) = integrate_history(modulus1_old.at(index),
                                                   sphere1_old[l].at(index),
                                                   sphere1_eval,
                                                   last_time_step_size);

            modulus2.at(index) = integrate_history(modulus2_old.at(index),
                                                   sphere2_old[l].at(index),
                                                   sphere2_eval,
                                                   last_time_step_size);

            modulus3.at(index) = integrate_history(modulus3_old.at(index),
                                                   sphere3_old[l].at(index),
                                                   sphere3_eval,
                                                   last_time_step_size);

            sphere1[l].at(index) = sphere1_eval;
            sphere2[l].at(index) = sphere2_eval;
            sphere3[l].at(index) = sphere3_eval;
        });

        matrix3 const macro_stress = compute_macro_stress(F_bar,
                                                          modulus1,
                                                          modulus2,
                                                          active_segments[l],
                                                          reduction_factor[l]);

        matrix6 const macro_moduli = compute_macro_moduli(F_bar,
                                                          modulus2,
                                                          modulus3,
                                                          active_segments[l],
                                                          reduction_factor[l]);

        auto const J = det_F[l];

        auto const pressure = J * volumetric_free_energy_dJ(J, K);

        cauchy_stresses[l] = compute_kirchhoff_stress(pressure, macro_stress) / J;

        tangent_operators[l] = compute_material_tangent(J, K, macro_moduli, macro_stress);
    });
}

matrix3 gaussian_ageing_affine_microsphere::compute_macro_stress(matrix3 const& F_bar,
                                                                 std::vector<double> const& modulus1,
                                                                 std::vector<double> const& modulus2,
                                                                 double const active_segments,
                                                                 double const reduction_factor) const
{
    return 3.0 * reduction_factor
           * unit_sphere
                 .integrate(matrix3::Zero().eval(), [&](auto const& coordinates, auto const index) -> matrix3 {
                     auto const& [r, _] = coordinates;

                     vector3 const t = deformed_tangent(F_bar, r);

                     auto const micro_stretch = compute_microstretch(t);

                     auto const shear_modulus = material.shear_modulus() + modulus1.at(index)
                                                - 0.5
                                                      * compute_prefactor_first(micro_stretch,
                                                                                active_segments)
                                                      * modulus2.at(index);

                     return shear_modulus * outer_product(t, t);
                 });
}

matrix6 gaussian_ageing_affine_microsphere::compute_macro_moduli(matrix3 const& F_bar,
                                                                 std::vector<double> const& modulus2,
                                                                 std::vector<double> const& modulus3,
                                                                 double const active_segments,
                                                                 double const reduction_factor) const
{
    return 3.0 * reduction_factor
           * unit_sphere
                 .integrate(matrix6::Zero().eval(), [&](auto const& coordinates, auto const index) -> matrix6 {
                     auto const& [r, _] = coordinates;

                     vector3 const t = deformed_tangent(F_bar, r);

                     auto const micro_stretch = compute_microstretch(t);

                     auto const shear_modulus = -1.0 / 2.0
                                                    * compute_prefactor_second(micro_stretch,
                                                                               active_segments)
                                                    * modulus2.at(index)
                                                + compute_prefactor_first(micro_stretch,
                                                                          active_segments)
                                                      * modulus3.at(index);

                     return shear_modulus * std::pow(micro_stretch, -1) * outer_product(t, t, t, t);
                 });
}

double gaussian_ageing_affine_microsphere::evaluate_integrand_first(double const creation_rate,
                                                                    double const reduction_factor,
                                                                    double const micro_stretch) const
{
    return creation_rate / (micro_stretch * reduction_factor);
}

double gaussian_ageing_affine_microsphere::evaluate_integrand_second(double const creation_rate,
                                                                     double const reduction_factor,
                                                                     double const micro_stretch) const
{
    return creation_rate / (std::pow(micro_stretch, 2) * reduction_factor);
}

double gaussian_ageing_affine_microsphere::evaluate_integrand_third(double const creation_rate,
                                                                    double const reduction_factor,
                                                                    double const micro_stretch) const
{
    return creation_rate / (std::pow(micro_stretch, 3) * reduction_factor);
}

double gaussian_ageing_affine_microsphere::compute_prefactor_first(double const micro_stretch,
                                                                   double const active_segments) const
{
    auto const b = material.bond_length();
    // auto const N = material.segments_per_chain();
    auto const N = active_segments;

    return micro_stretch
           - std::pow(micro_stretch, -1) * std::log(3.0 / (2.0 * M_PI * N * std::pow(b, 2)));
}

double gaussian_ageing_affine_microsphere::compute_prefactor_second(double const micro_stretch,
                                                                    double const active_segments) const
{
    auto const b = material.bond_length();
    // auto const N = material.segments_per_chain();
    auto const N = active_segments;

    return 3.0 + std::pow(micro_stretch, -2) * std::log(3.0 / (2.0 * M_PI * N * std::pow(b, 2)));
}

double gaussian_ageing_affine_microsphere::integrate_history(double const modulus,
                                                             double const sphere_old,
                                                             double const sphere,
                                                             double const time_step_size) const
{
    return modulus + partial_trapezoidal(sphere_old, sphere, time_step_size);
}
}
