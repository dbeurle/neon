
#include "constitutive/mechanical/solid/gaussian_ageing_affine_microsphere.hpp"

#include "constitutive/internal_variables.hpp"
#include "constitutive/mechanical/detail/microsphere.hpp"
#include "constitutive/mechanical/volumetric_free_energy.hpp"

#include <tbb/parallel_for.h>

#include <numeric>

using neon::mechanical::solid::gaussian_ageing_affine_microsphere;

gaussian_ageing_affine_microsphere::gaussian_ageing_affine_microsphere(
    std::shared_ptr<internal_variables_t>& variables,
    json const& material_data,
    unit_sphere_quadrature::Rule const rule)
    : gaussian_affine_microsphere(variables, material_data, rule), material(material_data)
{
    // Initialise the history variables for the chemical ageing network
    secondary_network.resize(variables->entries(), {1, matrix3::Identity()});

    // segments.resize(variables->entries(), {1, material.segments_per_chain()});

    cross_link_density.resize(variables->entries(), {1, material.shear_modulus()});
}

void gaussian_ageing_affine_microsphere::update_internal_variables(double const time_step_size)
{
    auto& tangent_operators = variables->fetch(internal_variables_t::rank4::tangent_operator);

    auto const& deformation_gradients = variables->fetch(
        internal_variables_t::Tensor::DeformationGradient);

    auto& cauchy_stresses = variables->fetch(internal_variables_t::Tensor::Cauchy);

    auto const& det_deformation_gradients = variables->fetch(internal_variables_t::scalar::DetF);

    auto const K{material.bulk_modulus()};
    auto const G{material.shear_modulus()};

    tbb::parallel_for(std::size_t{0}, deformation_gradients.size(), [&, this](auto const l) {
        secondary_network[l].push_back(deformation_gradients[l]);

        auto const& F_history = secondary_network[l];

        // Compute the deviatoric (unimodular) deformation gradient
        matrix3 const F_bar = unimodular(deformation_gradients[l]);

        auto const J = det_deformation_gradients[l];

        auto const pressure = J * volumetric_free_energy_dJ(J, K);

        // Project the stresses to obtain the Kirchhoff macro-stress

        // The stress is a computation over the deformation history
        matrix3 const macro_stress = std::accumulate(std::begin(F_history),
                                                     std::end(F_history),
                                                     matrix3::Zero().eval(),
                                                     [h = 0ul](matrix3 sum, auto const& F) mutable {
                                                         ++h;
                                                         return sum + compute_macro_stress(F_bar, G);
                                                     });

        matrix3 const macro_stress = compute_macro_stress(F_bar, G);

        cauchy_stresses[l] = compute_kirchhoff_stress(pressure, macro_stress) / J;

        tangent_operators[l] = compute_material_tangent(J,
                                                        K,
                                                        compute_macro_moduli(F_bar, G),
                                                        macro_stress);
    });
}