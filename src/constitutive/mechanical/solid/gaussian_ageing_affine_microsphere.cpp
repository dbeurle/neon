
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
    auto n{0};
    for (auto const& d : shear_moduli)
    {
        std::cout << "Quadrature point " << n << "\n";
        for (auto i : d)
        {
            std::cout << i << "\n";
        }
        ++n;
    }
}

void gaussian_ageing_affine_microsphere::update_internal_variables(double const time_step_size)
{
    auto& tangent_operators = variables->fetch(internal_variables_t::rank4::tangent_operator);

    auto const& deformation_gradients = variables->fetch(
        internal_variables_t::Tensor::DeformationGradient);

    auto& cauchy_stresses = variables->fetch(internal_variables_t::Tensor::Cauchy);

    auto const& det_deformation_gradients = variables->fetch(internal_variables_t::scalar::DetF);

    auto const K{material.bulk_modulus()};

    tbb::parallel_for(std::size_t{0}, deformation_gradients.size(), [&, this](auto const l) {
        matrix3 const& F = deformation_gradients[l];

        if (!is_approx(time_step_size, 0.0))
        {
            if (l == 0)
            {
                std::cout << "Shear moduli before\n";
                for (auto i : shear_moduli[l]) std::cout << i << "\n";
                std::cout << "\n";
            }

            // Update the network due to softening (chain scission)
            shear_moduli[l] = material.compute_shear_moduli(shear_moduli[l],
                                                            segments[l],
                                                            time_step_size);

            if (l == 0)
            {
                std::cout << "Shear moduli after\n";
                for (auto i : shear_moduli[l]) std::cout << i << "\n";
                std::cout << "\n";
            }

            // Add a new entry at the for the newly formed secondary network
            shear_moduli[l].emplace_back(material.compute_new_shear_modulus(time_step_size));

            if (l == 0)
            {
                std::cout << "Shear moduli update\n";
                for (auto i : shear_moduli[l]) std::cout << i << "\n";
                std::cout << "\n";
            }

            segments[l].emplace_back(material.compute_new_segment(segments[l].back(), time_step_size));

            if (l == 0)
            {
                std::cout << "Segments after update\n";
                for (auto i : segments[l]) std::cout << i << "\n";
                std::cout << "\n";
            }

            intermediate_deformations[l].push_back(F);

            if (l == 0)
            {
                std::cout << "Deformation gradient after update\n";
                for (auto i : intermediate_deformations[l]) std::cout << i << "\n";
                std::cout << "\n";
            }
        }

        // Stress computation over the secondary network deformation history
        // clang-format off
        matrix3 const macro_stress = std::accumulate(std::begin(intermediate_deformations[l]),
                                                     std::end(intermediate_deformations[l]),
                                                     matrix3::Zero().eval(),
                                                     [&, this, h = 0ul](matrix3 const p, auto const& F_0) mutable {

                                                         matrix3 const F_bar = unimodular(F_0.inverse() * F);

                                                         auto const G{shear_moduli[l][h++]};

                                                         return p + compute_macro_stress(F_bar, G);
                                                     });

        // Tangent operator computation over the secondary network deformation history
        matrix6 const macro_moduli = std::accumulate(std::begin(intermediate_deformations[l]),
                                                     std::end(intermediate_deformations[l]),
                                                     matrix6::Zero().eval(),
                                                     [&, this, h = 0ul](matrix6 const p, auto const& F_0) mutable {
                                                         matrix3 const F_bar = unimodular(F_0.inverse() * F);

                                                         auto const G{shear_moduli[l][h++]};

                                                         return p + compute_macro_moduli(F_bar, G);
                                                     });
        // clang-format on

        auto const J = det_deformation_gradients[l];

        auto const pressure = J * volumetric_free_energy_dJ(J, K);

        cauchy_stresses[l] = compute_kirchhoff_stress(pressure, macro_stress) / J;

        tangent_operators[l] = compute_material_tangent(J, K, macro_moduli, macro_stress);
    });
}
