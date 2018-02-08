
#include "gaussian_affine_microsphere_incremental.hpp"

#include "constitutive/InternalVariables.hpp"
#include "constitutive/mechanical/volumetric_free_energy.hpp"

#include "numeric/dense_matrix.hpp"

#include <range/v3/view/transform.hpp>
#include <range/v3/view/zip.hpp>

#include <omp.h>

namespace neon::mechanical::solid
{
void gaussian_affine_microsphere_incremental::update_internal_variables(double const time_step_size)
{
    using ranges::view::transform;
    using ranges::view::zip;

    auto& tangent_operators = variables->fetch(InternalVariables::rank4::tangent_operator);

    auto const& deformation_gradients = variables->fetch(
        InternalVariables::Tensor::DeformationGradient);

    auto const& deformation_gradients_old = variables->fetch_old(
        InternalVariables::Tensor::DeformationGradient);

    auto& cauchy_stresses = variables->fetch(InternalVariables::Tensor::Cauchy);
    auto& macro_stresses = variables->fetch(InternalVariables::Tensor::Kirchhoff);

    auto const& det_deformation_gradients = variables->fetch(InternalVariables::Scalar::DetF);

    auto const K = material.bulk_modulus();
    auto const G = material.shear_modulus();
    auto const N = material.segments_per_chain();

    // Compute the incremental deformation gradient
    auto const incremental_deformation_gradients = zip(deformation_gradients,
                                                       deformation_gradients_old)
                                                   | transform([&](auto const& tpl) -> matrix3 {
                                                         auto const & [ F, F_old ] = tpl;

                                                         return F * F_old.inverse();
                                                     });

// Compute the macrostresses on the unit sphere
#pragma omp parallel for
    for (auto l = 0; l < deformation_gradients.size(); ++l)
    {
        auto const& F_inc = incremental_deformation_gradients[l];
        macro_stresses[l] = compute_macro_stress(unimodular(F_inc), G, N);
    }

    // Project the stresses to obtain the Cauchy stress
    cauchy_stresses = zip(macro_stresses, det_deformation_gradients)
                      | transform([&](auto const& tpl) -> matrix3 {
                            auto const & [ macro_stress, J ] = tpl;

                            auto const pressure = J * volumetric_free_energy_dJ(J, K);

                            return compute_kirchhoff_stress(pressure, macro_stress) / J;
                        });

#pragma omp parallel for
    for (auto l = 0; l < deformation_gradients.size(); ++l)
    {
        auto const& F = deformation_gradients[l];
        auto const& macro_stress = macro_stresses[l];
        auto const& J = det_deformation_gradients[l];

        tangent_operators[l] = compute_material_tangent(J,
                                                        K,
                                                        compute_macro_moduli(unimodular(F), G, N),
                                                        macro_stress);
    }
}
}
