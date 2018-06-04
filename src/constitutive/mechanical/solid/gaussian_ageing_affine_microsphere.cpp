
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

using neon::mechanical::solid::gaussian_ageing_affine_microsphere;

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

    variables->add(internal_variables_t::vector::sphere_integral,
                   internal_variables_t::vector::sphere_last_function);

    auto [sphere_integral,
          sphere_last_function] = variables->get(internal_variables_t::vector::sphere_integral,
                                                 internal_variables_t::vector::sphere_last_function);

    for (auto& i : sphere_integral)
    {
        i.resize(unit_sphere.points(), 0.0);
    }
    for (auto& i : sphere_last_function)
    {
        i.resize(unit_sphere.points(), 0.0);
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
    std::cout << "Gaussian ageing model!" << std::endl;

    auto& tangent_operators = variables->get(internal_variables_t::fourth::tangent_operator);

    auto const& deformation_gradients = variables->get(
        internal_variables_t::second::DeformationGradient);

    auto& cauchy_stresses = variables->get(internal_variables_t::second::cauchy_stress);

    auto const& det_F = variables->get(internal_variables_t::scalar::DetF);

    auto [active_shear_modulus,
          inactive_shear_modulus,
          reduction_factor,
          active_segments,
          inactive_segments] = variables->get(internal_variables_t::scalar::active_shear_modulus,
                                              internal_variables_t::scalar::inactive_shear_modulus,
                                              internal_variables_t::scalar::reduction_factor,
                                              internal_variables_t::scalar::active_segments,
                                              internal_variables_t::scalar::inactive_segments);

    auto const& active_shear_modulus_old = variables->get_old(
        internal_variables_t::scalar::active_shear_modulus);

    auto [sphere_integral,
          sphere_last_function] = variables->get(internal_variable_t::vector::sphere_integral,
                                                 internal_variable_t::vector::sphere_last_function);

    auto const K{material.bulk_modulus()};

    tbb::parallel_for(std::size_t{0}, deformation_gradients.size(), [&, this](auto const l) {
        matrix3 const& F = deformation_gradients[l];

        matrix3 const F_bar = unimodular(F);

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

            // Partially integrate the ageing integral and accumulate into the
            // previous integral.  Save the current value for the next step in
            // the accumulated integral
            std::cout << "parameters\n" << parameters << std::endl;
            std::cout << "material.creation_rate(parameters) " << material.creation_rate(parameters)
                      << std::endl;

            // Loop over all the pointers
            unit_sphere.for_each([&](auto const& coordinates, auto const sphere_index) {
                auto const& [r, _] = coordinates;

                auto const micro_stretch = compute_microstretch(deformed_tangent(F_bar, r));

                auto const sphere_function = material.creation_rate(parameters)
                                             / (parameters(2) * micro_stretch);

                sphere_integral.at(l).at(
                    sphere_index) += partial_trapezoidal(sphere_last_function.at(l).at(sphere_index),
                                                         sphere_function,
                                                         time_step_size);
                sphere_last_function.at(l).at(sphere_index) = sphere_function;
            });
        }
        // clang-format off
        matrix3 const macro_stress = 3.0 * unit_sphere.integrate(matrix3::Zero().eval(), [&](auto const& coordinates, auto const sphere_index) -> matrix3 {
                                            auto const& [r, _] = coordinates;

                                            vector3 const t = deformed_tangent(F_bar, r);

                                            auto const secondary_shear_modulus = sphere_integral.at(l).at(sphere_index);

                                            auto const shear_modulus = (material.shear_modulus() + secondary_shear_modulus) * reduction_factor[l];

                                            return shear_modulus * outer_product(t, t);
                                        });



        matrix6 const macro_moduli = -3.0 * unit_sphere.integrate(matrix6::Zero().eval(), [&](auto const& coordinates, auto const sphere_index) -> matrix6 {
                                                auto const& [r, _] = coordinates;

                                                vector3 const t = deformed_tangent(F_bar, r);

                                                auto const micro_stretch = compute_microstretch(t);

                                                auto const secondary_shear_modulus = sphere_integral.at(l).at(sphere_index);

                                                auto const shear_modulus = (material.shear_modulus() + secondary_shear_modulus) * reduction_factor[l];

                                                return shear_modulus * std::pow(micro_stretch, -2) * outer_product(t, t, t, t);
                                            });
        // clang-format on

        auto const J = det_F[l];

        auto const pressure = J * volumetric_free_energy_dJ(J, K);

        cauchy_stresses[l] = compute_kirchhoff_stress(pressure, macro_stress) / J;

        tangent_operators[l] = compute_material_tangent(J, K, macro_moduli, macro_stress);
    });
}
