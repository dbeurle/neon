
#include "NonAffineMicrosphere.hpp"

#include "InternalVariables.hpp"
#include "numeric/DenseTypes.hpp"

#include <json/json.h>
#include <range/v3/algorithm.hpp>
#include <range/v3/view.hpp>

#include <exception>
#include <iostream>
#include <omp.h>

namespace neon
{
NonAffineMicrosphere::NonAffineMicrosphere(InternalVariables& variables,
                                           Json::Value const& material_data)
    : AffineMicrosphere(variables, material_data)
{
    if (!material_data.isMember("NonAffineStretchParameter"))
    {
        throw std::runtime_error("\"NonAffineStretchParameter\" not specified in material data\n");
    }
    non_affine_stretch_parameter = material_data["NonAffineStretchParameter"].asDouble();
}

void NonAffineMicrosphere::update_internal_variables(double const time_step_size)
{
    using ranges::view::transform;
    using ranges::view::zip;

    // TODO Change these back once OpenMP allows structured bindings

    // Get references into the hash table
    auto& F_list = variables(InternalVariables::Tensor::DeformationGradient);
    auto& cauchy_stress_list = variables(InternalVariables::Tensor::Cauchy);
    auto& dev_stress_list = variables(InternalVariables::Tensor::Kirchhoff);

    auto const& detF_list = variables(InternalVariables::Scalar::DetF);

    // Material properties
    auto const K = material.bulk_modulus();
    auto const G_eff = material.shear_modulus();
    auto const N = material.segments_per_chain();

    auto const p = non_affine_stretch_parameter;

#pragma omp parallel for
    for (auto l = 0; l < F_list.size(); ++l)
    {
        auto& stress_dev = dev_stress_list[l]; // Deviatoric stress
        auto const& F = F_list[l];             // Deformation gradient
        auto const& J = detF_list[l];          // Determinant of the deformation gradient

        Matrix3 const F_deviatoric = std::pow(J, -1.0 / 3.0) * F;

        // Compute the non-affine stretch parts
        auto const nonaffine_stretch = compute_nonaffine_stretch(F_deviatoric);

        Matrix3 const h = compute_h_tensor(F_deviatoric);
        Matrix6 const H = compute_H_tensor(F_deviatoric);

        // Compute the microstress and micro moduli
        auto const micro_kirchhoff_f = G_eff * (3.0 * N - std::pow(nonaffine_stretch, 2))
                                       / ((N - std::pow(nonaffine_stretch, 2)) * nonaffine_stretch);

        auto const micro_moduli_f = G_eff * (std::pow(nonaffine_stretch, 4) + 3.0 * std::pow(N, 2))
                                    / std::pow(N - std::pow(nonaffine_stretch, 2), 2);

        // Compute the macrostress and macromoduli for chain force
        Matrix3 const macro_kirchhoff_f = micro_kirchhoff_f * std::pow(nonaffine_stretch, 1.0 - p)
                                          * h;

        Matrix6 const macro_moduli_f = (micro_moduli_f * std::pow(nonaffine_stretch, 2.0 - 2.0 * p)
                                        - (p - 1.0) * micro_kirchhoff_f
                                              * std::pow(nonaffine_stretch, 1.0 - 2.0 * p))
                                           * outer_product(h, h)
                                       + micro_kirchhoff_f * std::pow(nonaffine_stretch, 1.0 - p) * H;

        // Compute the non-affine tube contribution
        Matrix3 const k = compute_k_tensor(F_deviatoric);

        Matrix6 const K = compute_K_tensor(F_deviatoric);
        Matrix6 const G = compute_G_tensor(F_deviatoric);

        // Compute the macrostress and macromoduli for tube contraint
        Matrix3 const macro_kirchhoff_c = -G_eff * N * effective_tube_geometry * k;

        Matrix6 const macro_moduli_c = G_eff * N * effective_tube_geometry * (K + G);

        // Superimposed stress response from tube and chain contributions
        Matrix3 const macro_kirchhoff = macro_kirchhoff_f + macro_kirchhoff_c;

        Matrix6 const macro_moduli = macro_moduli_f + macro_moduli_c;

        // Perform the deviatoric projection for the stress and moduli
    }
}

double NonAffineMicrosphere::compute_nonaffine_stretch(Matrix3 const& F_deviatoric) const
{
    auto const p = non_affine_stretch_parameter;

    return std::pow(unit_sphere.integrate(0.0,
                                          [&](auto const& coordinates, auto const& l) {
                                              auto const & [ r, r_outer_r ] = coordinates;

                                              // Deformed tangent
                                              Vector3 const t = F_deviatoric * r;

                                              return std::pow(t.norm(), p);
                                          }),
                    1.0 / p);
}

Matrix3 NonAffineMicrosphere::compute_h_tensor(Matrix3 const& F_deviatoric) const
{
    auto const p = non_affine_stretch_parameter;

    return unit_sphere.integrate(Matrix3::Zero().eval(), [&](auto const& xyz, auto const& l) {
        auto const & [ r, r_o_r ] = xyz;

        Vector3 const t = F_deviatoric * r;

        return std::pow(t.norm(), p - 2.0) * outer_product(t, t);
    });
}

Matrix6 NonAffineMicrosphere::compute_H_tensor(Matrix3 const& F_deviatoric) const
{
    auto const p = non_affine_stretch_parameter;

    return (p - 2.0)
           * unit_sphere.integrate(Matrix6::Zero().eval(), [&](auto const& xyz, auto const& l) {
                 auto const & [ r, r_o_r ] = xyz;

                 Vector3 const t = F_deviatoric * r;

                 return std::pow(t.norm(), p - 4.0) * outer_product(t, t, t, t);
             });
}

Matrix3 NonAffineMicrosphere::compute_k_tensor(Matrix3 const& F_deviatoric) const
{
    auto const q = non_affine_tube_parameter;

    return unit_sphere.integrate(Matrix3::Zero().eval(), [&](auto const& xyz, auto const& l) {
        auto const & [ r, r_o_r ] = xyz;

        // Compute deformed normal
        Vector3 const n = F_deviatoric.inverse().transpose() * r;

        return std::pow(n.norm(), q - 2.0) * outer_product(n, n);
    });
}

Matrix6 NonAffineMicrosphere::compute_K_tensor(Matrix3 const& F_deviatoric) const
{
    auto const q = non_affine_tube_parameter;

    return q * (q - 2.0)
           * unit_sphere.integrate(Matrix6::Zero().eval(), [&](auto const& xyz, auto const& l) {
                 auto const & [ r, r_o_r ] = xyz;

                 // Compute deformed normal
                 Vector3 const n = F_deviatoric.inverse().transpose() * r;

                 return std::pow(n.norm(), q - 4.0) * outer_product(n, n, n, n);
             });
}

Matrix6 NonAffineMicrosphere::compute_G_tensor(Matrix3 const& F_deviatoric) const
{
    auto const q = non_affine_tube_parameter;

    return 2.0 * q
           * unit_sphere.integrate(Matrix6::Zero().eval(), [&](auto const& xyz, auto const& l) {
                 auto const & [ r, r_o_r ] = xyz;

                 // Compute deformed normal
                 Vector3 const n = F_deviatoric.inverse().transpose() * r;

                 return std::pow(n.norm(), q - 2.0) * compute_o_dot_product(n);
             });
}
}
