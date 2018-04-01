
#include "constitutive/mechanical/solid/gaussian_ageing_affine_microsphere.hpp"

#include "constitutive/internal_variables.hpp"
#include "constitutive/mechanical/detail/microsphere.hpp"
#include "constitutive/mechanical/volumetric_free_energy.hpp"
#include "solver/time/runge_kutta_integration.hpp"
#include "numeric/float_compare.hpp"

#include <tbb/parallel_for.h>

#include <cassert>
#include <numeric>
#include <iostream>

using neon::mechanical::solid::gaussian_ageing_affine_microsphere;

gaussian_ageing_affine_microsphere::gaussian_ageing_affine_microsphere(
    std::shared_ptr<internal_variables_t>& variables,
    json const& material_data,
    unit_sphere_quadrature::Rule const rule)
    : gaussian_affine_microsphere{variables, material_data, rule},
      material{material_data},
      shear_moduli{variables->entries(), {material.shear_modulus()}},
      inactive_shear_moduli(variables->entries(), 0.0),
      segments{variables->entries(), {material.segments_per_chain(), 0.0}},
      intermediate_deformations{variables->entries(), {matrix3::Identity()}}
{
    variables->add(internal_variables_t::scalar::active_chains);
    variables->add(internal_variables_t::scalar::inactive_chains);

    variables->add(internal_variables_t::scalar::active_segment_average);
    variables->add(internal_variables_t::scalar::inactive_segment_average);

    auto& active_chains = variables->fetch(internal_variables_t::scalar::active_chains);

    std::fill(begin(active_chains), end(active_chains), material.shear_modulus());

    variables->commit();
}

void gaussian_ageing_affine_microsphere::update_internal_variables(double const time_step_size)
{
    auto& tangent_operators = variables->fetch(internal_variables_t::rank4::tangent_operator);

    auto const& deformation_gradients = variables->fetch(
        internal_variables_t::Tensor::DeformationGradient);

    auto& cauchy_stresses = variables->fetch(internal_variables_t::Tensor::Cauchy);

    auto const& det_deformation_gradients = variables->fetch(internal_variables_t::scalar::DetF);

    auto& active_chains = variables->fetch(internal_variables_t::scalar::active_chains);

    auto const K{material.bulk_modulus()};

    tbb::parallel_for(std::size_t{0}, deformation_gradients.size(), [&, this](auto const l) {
        matrix3 const& F = deformation_gradients[l];

        if (!is_approx(time_step_size, 0.0))
        {
            // Perform the updates
            auto const p_s = material.scission_probability();
            auto const p_c = material.recombination_probability();

            auto& shear_moduli_history = shear_moduli[l];

            // Create a new secondary network
            shear_moduli_history.push_back(0.0);

            // Unpack the segments
            auto& [N_a, N_ia] = segments[l];

            // Define the right hand side of the ageing evolution equations
            auto integrator = runge_kutta_fourth_order([&](auto const t, vector y) -> vector {
                assert(y.size() >= 5);

                auto const active_set_count = y.size() - 3;

                // Extract the previous active sets
                vector const n_a = y.head(active_set_count - 1);

                // Inactive set
                double const n_ia = y(active_set_count);

                // Average number of segments in active set
                auto const N_a = y(active_set_count + 1);

                // Average number of segments in inactive set
                auto const N_ia = y(active_set_count + 2);

                // Inactive set recombination
                auto const alpha = 1.0 - std::pow(1.0 - p_c, N_ia + 1.0);
                // Active set scission
                auto const beta = 1.0 - std::pow(1.0 - p_s, N_a);
                // Active set generation
                auto const eta = 1.0 - std::pow(1.0 - p_c, N_a + 1.0);
                // Inactive set generation
                auto const nu = 1.0 - std::pow(1.0 - p_s, N_ia);

                vector const f1 = -n_a * (beta + 2.0 * eta);

                // Latest active set created
                auto const f2 = alpha * n_ia + 4.0 * eta * n_a.sum();

                // Inactive set rate of change
                auto const f3 = 2.0 * (beta * n_a.sum() + n_ia * (nu - alpha));

                // Total rate of change for the active sets (d (n_a) / dt)
                auto const f4 = (2.0 * eta - beta) * n_a.sum() + alpha * n_ia;

                // Active set average segments per chain rate
                auto const g1 = (-N_a * beta * n_a.sum() + 2.0 * alpha * n_ia * N_ia - N_a * f4)
                                / n_a.sum();

                // Inactive set average segments per chain rate
                auto const g2 = n_ia > 1.0e-8
                                    ? (beta * N_a * n_a.sum() - 2.0 * alpha * N_ia * n_ia - N_ia * f3)
                                          / n_ia
                                    : 0.0;

                // Repack the data into y and return this as the update
                y.head(f1.size()) = f1;
                y(f1.size()) = f2;
                // Inactive set
                y(f1.size() + 1) = f3;
                // Segments
                y(f1.size() + 2) = g1;
                y(f1.size() + 3) = g2;

                return y;
            });

            // Fill a Eigen vector (active sets, inactive set, segments (active and inactive))
            vector z(shear_moduli_history.size() + 3);

            std::copy(begin(shear_moduli_history), end(shear_moduli_history), z.data());

            z(z.size() - 3) = inactive_shear_moduli[l];
            z(z.size() - 2) = N_a;
            z(z.size() - 1) = N_ia;

            // Integrate the system
            z += integrator(0.0, z, time_step_size).eval();

            for (std::int64_t i{}; i < z.size(); ++i)
            {
                if (z(i) < 0.0) z(i) = 0.0;
            }

            std::transform(begin(shear_moduli_history),
                           end(shear_moduli_history),
                           begin(shear_moduli_history),
                           [&z, index = 0ul](auto const i) mutable { return z(index++); });

            // Use z to update the results
            inactive_shear_moduli[l] = z(z.size() - 3);
            N_a = z(z.size() - 2);
            N_ia = z(z.size() - 1);

            intermediate_deformations[l].push_back(F);

            // Update the history variables for plotting
            variables->fetch(internal_variables_t::scalar::active_chains)[l] = std::
                accumulate(begin(shear_moduli_history), end(shear_moduli_history), 0.0);
            variables->fetch(internal_variables_t::scalar::inactive_chains)[l] = z(z.size() - 3);
            variables->fetch(internal_variables_t::scalar::active_segment_average)[l] = N_a;
            variables->fetch(internal_variables_t::scalar::inactive_segment_average)[l] = N_ia;
        }
        // Stress computation over the secondary network deformation history
        // clang-format off
        matrix3 const macro_stress = std::accumulate(begin(intermediate_deformations[l]),
                                                     end(intermediate_deformations[l]),
                                                     matrix3::Zero().eval(),
                                                     [&, this, time_index = 0ul](matrix3 const p, auto const& F_0) mutable {

                                                         matrix3 const F_bar = unimodular(F_0.inverse() * F);

                                                         auto const G{shear_moduli[l][time_index++]};

                                                         return p + compute_macro_stress(F_bar, G);
                                                     });

        // Tangent operator computation over the secondary network deformation history
        matrix6 const macro_moduli = std::accumulate(begin(intermediate_deformations[l]),
                                                     end(intermediate_deformations[l]),
                                                     matrix6::Zero().eval(),
                                                     [&, this, time_index = 0ul](matrix6 const p, auto const& F_0) mutable {

                                                         matrix3 const F_bar = unimodular(F_0.inverse() * F);

                                                         auto const G{shear_moduli[l][time_index++]};

                                                         return p + compute_macro_moduli(F_bar, G);
                                                     });
        // clang-format on
        auto const J = det_deformation_gradients[l];

        auto const pressure = J * volumetric_free_energy_dJ(J, K);

        cauchy_stresses[l] = compute_kirchhoff_stress(pressure, macro_stress) / J;

        tangent_operators[l] = compute_material_tangent(J, K, macro_moduli, macro_stress);
    });
}
