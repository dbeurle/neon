
#include "constitutive/mechanics/solid/gaussian_ageing_affine_microsphere.hpp"

#include "constitutive/internal_variables.hpp"
#include "constitutive/mechanics/detail/microsphere.hpp"
#include "constitutive/mechanics/volumetric_free_energy.hpp"
#include "quadrature/trapezoidal.hpp"
#include "solver/time/runge_kutta_integration.hpp"
#include "numeric/float_compare.hpp"

#include <tbb/parallel_for.h>

#include <cassert>
#include <numeric>

namespace neon::mechanics::solid
{
gaussian_ageing_affine_microsphere::gaussian_ageing_affine_microsphere(
    std::shared_ptr<internal_variables_t>& variables,
    json const& material_data,
    unit_sphere_quadrature::point const rule)
    : gaussian_affine_microsphere{variables, material_data, rule}, material{material_data}
{
    variables->add(variable::scalar::active_shear_modulus,
                   variable::scalar::inactive_shear_modulus,
                   variable::scalar::active_segments,
                   variable::scalar::inactive_segments,
                   variable::scalar::reduction_factor);

    variables->add(variable::vector::accumulated_ageing_integral,
                   variable::vector::previous_integrand);

    names.emplace("active_shear_modulus");
    names.emplace("inactive_shear_modulus");
    names.emplace("active_segments");
    names.emplace("inactive_segments");
    names.emplace("reduction_factor");

    for (auto& values : variables->get(variable::vector::accumulated_ageing_integral))
    {
        values.resize(unit_sphere.points(), 0.0);
    }
    for (auto& values : variables->get(variable::vector::previous_integrand))
    {
        values.resize(unit_sphere.points(), 0.0);
    }

    auto [active_shear_modulus,
          active_segments,
          reduction] = variables->get(variable::scalar::active_shear_modulus,
                                      variable::scalar::active_segments,
                                      variable::scalar::reduction_factor);

    std::fill(begin(active_shear_modulus), end(active_shear_modulus), material.shear_modulus());
    std::fill(begin(active_segments), end(active_segments), material.segments_per_chain());
    std::fill(begin(reduction), end(reduction), 1.0);

    variables->commit();
}

void gaussian_ageing_affine_microsphere::update_internal_variables(double const time_step_size)
{
    // Polymer network values
    auto& active_shear_modulus = variables->get(variable::scalar::active_shear_modulus);
    auto& inactive_shear_modulus = variables->get(variable::scalar::inactive_shear_modulus);
    auto& reduction_factor = variables->get(variable::scalar::reduction_factor);
    auto& active_segments = variables->get(variable::scalar::active_segments);
    auto& inactive_segments = variables->get(variable::scalar::inactive_segments);

    // accumulated integral values from time integration
    auto& moduli = variables->get(variable::vector::accumulated_ageing_integral);
    auto& last_evaluation = variables->get(variable::vector::previous_integrand);

    auto const& moduli_old = variables->get_old(variable::vector::accumulated_ageing_integral);
    auto const& last_evaluation_old = variables->get_old(variable::vector::previous_integrand);

    auto const& deformation_gradients = variables->get(variable::second::deformation_gradient);
    auto const& det_F = variables->get(variable::scalar::DetF);

    auto& cauchy_stresses = variables->get(variable::second::cauchy_stress);
    auto& tangent_operators = variables->get(variable::fourth::tangent_operator);

    auto const K{material.bulk_modulus()};

    tbb::parallel_for(std::size_t{0}, deformation_gradients.size(), [&, this](auto const l) {
        // unimodular deformation gradient
        matrix3 const F_bar = unimodular(deformation_gradients[l]);

        auto& modulus = moduli[l];
        auto& last_h = last_evaluation[l];

        auto const& modulus_old = moduli_old[l];
        auto const& last_h_old = last_evaluation_old[l];

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

        auto const creation_rate = material.creation_rate(active_shear_modulus[l],
                                                          inactive_shear_modulus[l],
                                                          active_segments[l],
                                                          inactive_segments[l]);

        // Partially integrate secondary modulus
        unit_sphere.for_each([&](auto const& coordinates, auto const index) {
            auto const& [r, _] = coordinates;

            auto const stretch = compute_microstretch(deformed_tangent(F_bar, r));

            last_h[index] = evaluate_integrand(creation_rate, reduction_factor[l], stretch);

            modulus[index] = integrate_history(modulus_old[index],
                                               last_h_old[index],
                                               last_h[index],
                                               last_time_step_size);
        });

        matrix3 const macro_stress = compute_macro_stress(F_bar, modulus, reduction_factor[l]);

        matrix6 const macro_moduli = compute_macro_moduli(F_bar, creation_rate, last_time_step_size);

        auto const J = det_F[l];

        auto const pressure = J * volumetric_free_energy_dJ(J, K);

        cauchy_stresses[l] = compute_kirchhoff_stress(pressure, macro_stress) / J;

        tangent_operators[l] = compute_material_tangent(J, K, macro_moduli, macro_stress);
    });
}

matrix3 gaussian_ageing_affine_microsphere::compute_macro_stress(matrix3 const& F_bar,
                                                                 std::vector<double> const& modulus,
                                                                 double const reduction_factor) const
{
    return 3.0 * reduction_factor
           * unit_sphere.integrate(matrix3::Zero().eval(),
                                   [&](auto const& coordinates, auto const index) -> matrix3 {
                                       auto const& [r, _] = coordinates;

                                       vector3 const t = deformed_tangent(F_bar, r);

                                       return (material.shear_modulus() + modulus[index])
                                              * outer_product(t, t);
                                   });
}

matrix6 gaussian_ageing_affine_microsphere::compute_macro_moduli(matrix3 const& F_bar,
                                                                 double const creation_rate,
                                                                 double const time_step_size) const
{
    return -3.0 / 2.0 * creation_rate * time_step_size
           * unit_sphere.integrate(matrix6::Zero().eval(),
                                   [&](auto const& coordinates, auto) -> matrix6 {
                                       auto const& [r, _] = coordinates;

                                       vector3 const t = deformed_tangent(F_bar, r);

                                       auto const stretch = compute_microstretch(t);

                                       return std::pow(stretch, -3) * outer_product(t, t, t, t);
                                   });
}

double gaussian_ageing_affine_microsphere::evaluate_integrand(double const creation_rate,
                                                              double const reduction_factor,
                                                              double const micro_stretch) const
{
    return creation_rate / (micro_stretch * reduction_factor);
}

double gaussian_ageing_affine_microsphere::integrate_history(double const modulus,
                                                             double const sphere_old,
                                                             double const sphere,
                                                             double const time_step_size) const
{
    return modulus + partial_trapezoidal(sphere_old, sphere, time_step_size);
}
}
