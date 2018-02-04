
#include "affine_microsphere.hpp"

#include "constitutive/InternalVariables.hpp"
#include "numeric/dense_matrix.hpp"

#include <range/v3/view/transform.hpp>
#include <range/v3/view/zip.hpp>

#include <omp.h>

namespace neon::mechanical::solid
{
affine_microsphere::affine_microsphere(std::shared_ptr<InternalVariables>& variables,
                                     json const& material_data,
                                     unit_sphere_quadrature::Rule const rule)
    : ConstitutiveModel(variables), unit_sphere(rule), material(material_data)
{
    variables->add(InternalVariables::rank4::tangent_operator);

    // Deviatoric stress
    variables->add(InternalVariables::Tensor::Kirchhoff);

    // Commit these to history in case of failure on first time step
    variables->commit();
}

void affine_microsphere::update_internal_variables(double const time_step_size)
{
    auto& tangent_operators = variables->fetch(InternalVariables::rank4::tangent_operator);

    auto const& deformation_gradients = variables->fetch(
        InternalVariables::Tensor::DeformationGradient);
    auto& cauchy_stresses = variables->fetch(InternalVariables::Tensor::Cauchy);
    auto& macro_stresses = variables->fetch(InternalVariables::Tensor::Kirchhoff);

    auto const& det_deformation_gradients = variables->fetch(InternalVariables::Scalar::DetF);

    auto const K = material.bulk_modulus();
    auto const G = material.shear_modulus();
    auto const N = material.segments_per_chain();

// Compute the macrostresses on the unit sphere
#pragma omp parallel for
    for (auto l = 0; l < deformation_gradients.size(); ++l)
    {
        auto const& F = deformation_gradients[l];
        macro_stresses[l] = compute_macro_stress(unimodular(F), G, N);
    }

    // Project the stresses to obtain the Cauchy stress
    cauchy_stresses = ranges::view::zip(macro_stresses, det_deformation_gradients)
                      | ranges::view::transform([&](auto const& tpl) -> matrix3 {
                            auto const& [macro_stress, J] = tpl;

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

matrix3 affine_microsphere::compute_kirchhoff_stress(double const pressure,
                                                    matrix3 const& macro_stress) const
{
    // clang-format off
    return pressure * matrix3::Identity() + voigt::kinetic::from(P * voigt::kinetic::to(macro_stress));
    // clang-format on
}

matrix6 affine_microsphere::compute_material_tangent(double const J,
                                                    double const K,
                                                    matrix6 const& macro_C,
                                                    matrix3 const& macro_stress) const
{
    auto const pressure = J * volumetric_free_energy_dJ(J, K);
    auto const kappa = std::pow(J, 2) * volumetric_free_energy_second_d2J(J, K);

    // clang-format off
    matrix6 const D = macro_C
                    + 2.0 / 3.0 * macro_stress.trace() * voigt::kinematic::identity()
                    - 2.0 / 3.0 * (outer_product(macro_stress, matrix3::Identity()) +
                                   outer_product(matrix3::Identity(), macro_stress));

    // clang-format on
    return (kappa + pressure) * IoI - 2.0 * pressure * I + P * D * P;
}

matrix3 affine_microsphere::compute_macro_stress(matrix3 const& F_unimodular,
                                                double const bulk_modulus,
                                                double const N) const
{
    return bulk_modulus
           * unit_sphere.integrate(matrix3::Zero().eval(),
                                   [&](auto const& coordinates, auto const& l) -> matrix3 {
                                       auto const& [r, r_outer_r] = coordinates;

                                       vector3 const t = deformed_tangent(F_unimodular, r);

                                       return pade_first(compute_microstretch(t), N)
                                              * outer_product(t, t);
                                   });
}

matrix6 affine_microsphere::compute_macro_moduli(matrix3 const& F_unimodular,
                                                double const bulk_modulus,
                                                double const N) const
{
    // clang-format off
    return bulk_modulus * unit_sphere.integrate(matrix6::Zero().eval(),
                                                [&](auto const& coordinates, auto const& l) -> matrix6 {
                                                    auto const & [ r, r_outer_r ] = coordinates;

                                                    vector3 const t = deformed_tangent(F_unimodular, r);

                                                    auto const micro_stretch = compute_microstretch(t);

                                                    auto const a = std::pow(micro_stretch, -2) * (pade_second(micro_stretch, N) - pade_first(micro_stretch, N));

                                                    return a * outer_product(t, t, t, t);
                                                });
    // clang-format on
}

AffineMicrosphereWithDegradation::AffineMicrosphereWithDegradation(
    std::shared_ptr<InternalVariables>& variables,
    json const& material_data,
    unit_sphere_quadrature::Rule const rule)
    : affine_microsphere(variables, material_data, rule), material(material_data)
{
    variables->add(InternalVariables::Scalar::Chains, InternalVariables::Scalar::ShearModuli);

    // Shrink these down to the correct size
    variables->fetch(InternalVariables::Scalar::Chains).resize(material.groups(), 0.0);
    variables->fetch(InternalVariables::Scalar::ShearModuli).resize(material.groups(), 0.0);
    variables->fetch(InternalVariables::Scalar::Chains).shrink_to_fit();
    variables->fetch(InternalVariables::Scalar::ShearModuli).shrink_to_fit();

    // Fill the data with material properties using the material class
    variables->fetch(InternalVariables::Scalar::Chains) = material.chain_groups();
    variables->fetch(InternalVariables::Scalar::ShearModuli) = material.shear_moduli_groups();

    // Commit these to history in case of failure on first time step
    variables->commit();
}

void AffineMicrosphereWithDegradation::update_internal_variables(double const time_step_size)
{
    //     using ranges::view::transform;
    //     using ranges::view::zip;
    //
    //     auto& tangent_operators = variables(InternalVariables::rank4::tangent_operator);
    //
    //     auto const& deformation_gradients =
    //     variables(InternalVariables::Tensor::DeformationGradient); auto& cauchy_stresses =
    //     variables(InternalVariables::Tensor::Cauchy); auto& macro_stresses =
    //     variables(InternalVariables::Tensor::Kirchhoff);
    //
    //     auto const& det_deformation_gradients = variables(InternalVariables::Scalar::DetF);
    //
    //     auto& n_list = variables(InternalVariables::Scalar::Chains);
    //     auto& G_list = variables(InternalVariables::Scalar::ShearModuli);
    //
    //     // Update the material properties
    //     n_list = material.update_chains(n_list, time_step_size);
    //     G_list = material.compute_shear_moduli(n_list);
    //
    //     auto const K = material.bulk_modulus();
    //
    // // Compute the macrostresses on the unit sphere
    // #pragma omp parallel for
    //     for (auto l = 0; l < deformation_gradients.size(); ++l)
    //     {
    //         auto const& F = deformation_gradients[l]; // Deformation gradient
    //
    //         macro_stresses[l] = weighting(zip(G_list, , matrix3::Zero().eval(), [&](auto const&
    //         N) -> matrix3 {
    //             return compute_macro_stress(unimodular(F), N);
    //         });
    //     }
    //
    //     // Project the stresses to obtain the Cauchy stress
    //     cauchy_stresses = zip(macro_stresses, det_deformation_gradients)
    //                       | transform([&](auto const& tpl) -> matrix3 {
    //                             auto const & [ macro_stress, J ] = tpl;
    //
    //                             auto const pressure = J * volumetric_free_energy_dJ(J, K);
    //
    //                             return compute_kirchhoff_stress(pressure, macro_stress) / J;
    //                         });
    //
    // #pragma omp parallel for
    //     for (auto l = 0; l < deformation_gradients.size(); ++l)
    //     {
    //         auto const& F = deformation_gradients[l];
    //         auto const& macro_stress = macro_stresses[l];
    //         auto const& J = det_deformation_gradients[l];
    //
    //         matrix6 const macro_C = weighting(G_list, matrix6::Zero().eval(), [&](auto const& N)
    //         -> matrix6 {
    //             return compute_macro_moduli(unimodular(F), N);
    //         });
    //
    //         tangent_operators[l] = compute_material_tangent(J, K, macro_C, macro_stress);
    //     }
}
}
