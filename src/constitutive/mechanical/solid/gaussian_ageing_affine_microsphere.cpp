
#include "constitutive/mechanical/solid/gaussian_ageing_affine_microsphere.hpp"

#include "constitutive/internal_variables.hpp"
#include "constitutive/mechanical/detail/microsphere.hpp"
#include "constitutive/mechanical/volumetric_free_energy.hpp"
#include "numeric/float_compare.hpp"

#include <tbb/parallel_for.h>

#include <numeric>
#include <iostream>

using neon::mechanical::solid::gaussian_ageing_affine_microsphere;

gaussian_ageing_affine_microsphere::gaussian_ageing_affine_microsphere(
    std::shared_ptr<internal_variables_t>& variables,
    json const& material_data,
    unit_sphere_quadrature::Rule const rule)
    : gaussian_affine_microsphere{variables, material_data, rule},
      material{material_data},
      segments{variables->entries(), {material.segments_per_chain()}},
      shear_moduli{variables->entries(), {material.shear_modulus()}},
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

            auto const k = material.active_segment_decay_rate();

            auto& shear_moduli_pair_history = shear_moduli[l];
            auto& segment_pair_history = segments[l];

            // Define the right hand side of the ageing evolution equations
            auto integrator = runge_kutta_fourth_fifth_order([&](auto const t, auto const y) {
                auto const n_a = y(0);  // y(0) active number of chain
                auto const n_ia = y(1); // y(1) inactive number of chains
                auto const N_a = y(2);  // y(2) average number of segments in active set
                auto const N_ia = y(3); // y(3) average number of segments in inactive set

                // Probabilities of network events
                auto const alpha = 1.0 - std::pow(1.0 - p_c, N_ia + 1); // Inactive set recombination
                auto const beta = 1.0 - std::pow(1.0 - p_f, N_a);       // Active set scission
                auto const eta = 1.0 - std::pow(1.0 - p_c, N_a + 1);    // Active set generation
                auto const nu = 1.0 - std::pow(1.0 - p_f, N_ia);        // Inactive generation

                // Active set f(t, y)
                y(0) = 2 * (alpha * n_ia - beta * n_a + eta * n_a);

                // Inactive set f(t, y)
                y(1) = 2 * (beta * n_a - alpha * n_ia + nu * n_ia);

                // Active set average number of segments
                y(2) = -k * N_a;

                // Active set average number of segments
                y(3) = n_ia > 1.0e-8 ? -y(2) * n_a / n_ia + 2 * alpha * (N_ia - N_a)
                                           + 2 * beta * (N_a - N_ia) * n_a / n_ia
                                           - 2 * eta * N_a * n_a / n_ia - 2 * nu * N_ia
                                     : 0.0;

                return y;
            });

            // Create the new one
            // Update the network due to softening (chain scission)
            shear_moduli[l] = material.scission(shear_moduli[l], segments[l], time_step_size);

            // Add a new entry at the for the newly formed secondary network
            shear_moduli[l].emplace_back(material.compute_new_shear_modulus(time_step_size));

            std::cout << "Shear modulus: " << shear_moduli[l].back() << "\n";

            active_chains[l] = shear_moduli[l].back();

            segments[l].emplace_back(material.compute_new_segment(segments[l].back(), time_step_size));

            // Update the previously formed secondary network
            for (std::size_t i{}; i < shear_moduli_pair_history.size(); ++i)
            {
                auto const [N_a, N_ia] = segment_pair_history[i];
                auto const [n_a, n_ia] = shear_moduli_pair_history[i];

                vector4 y(4);
                y << n_a, n_ia, N_a, N_ia;

                // Integrate the system
                y += integrator(0.0, y, time_step_size);

                shear_moduli_pair_history[i] = {y(0), y(1)};
                segment_pair_history[i] = {y(2), y(3)};
            }

            intermediate_deformations[l].push_back(F);
        }
        // Stress computation over the secondary network deformation history
        // clang-format off
        matrix3 const macro_stress = std::accumulate(std::begin(intermediate_deformations[l]),
                                                     std::end(intermediate_deformations[l]),
                                                     matrix3::Zero().eval(),
                                                     [&, this, time_index = 0ul](matrix3 const p, auto const& F_0) mutable {

                                                         matrix3 const F_bar = unimodular(F_0.inverse() * F);

                                                         auto const G{shear_moduli[l][time_index++]};

                                                         return p + compute_macro_stress(F_bar, G);
                                                     });

        // Tangent operator computation over the secondary network deformation history
        matrix6 const macro_moduli = std::accumulate(std::begin(intermediate_deformations[l]),
                                                     std::end(intermediate_deformations[l]),
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
