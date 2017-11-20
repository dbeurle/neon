
#include "NonAffineMicrosphere.hpp"

#include "InternalVariables.hpp"
#include "numeric/DenseMatrix.hpp"

#include <json/value.h>

#include <exception>
#include <omp.h>

namespace neon::mech::solid
{
NonAffineMicrosphere::NonAffineMicrosphere(InternalVariables& variables,
                                           Json::Value const& material_data,
                                           UnitSphereQuadrature::Rule const rule)
    : AffineMicrosphere(variables, material_data, rule), material(material_data)
{
    if (!material_data.isMember("NonAffineStretchParameter"))
    {
        throw std::runtime_error("\"NonAffineStretchParameter\" not specified in material data\n");
    }
    non_affine_stretch_parameter = material_data["NonAffineStretchParameter"].asDouble();
}

void NonAffineMicrosphere::update_internal_variables(double const time_step_size)
{
    auto const& deformation_gradients = variables(InternalVariables::Tensor::DeformationGradient);
    auto& cauchy_stresses = variables(InternalVariables::Tensor::Cauchy);

    auto const& detF_list = variables(InternalVariables::Scalar::DetF);

    // Compute tangent moduli
    auto& tangent_operators = variables(InternalVariables::Matrix::TangentOperator);

    // Material properties
    auto const K_eff = material.bulk_modulus();
    auto const G_eff = material.shear_modulus();
    auto const N = material.segments_per_chain();
    auto const p = non_affine_stretch_parameter;

    auto const number_of_internal_variables = deformation_gradients.size();

#pragma omp parallel for
    for (auto l = 0; l < number_of_internal_variables; ++l)
    {
        auto const& F = deformation_gradients[l]; // Deformation gradient
        auto const& J = detF_list[l];             // Determinant of the deformation gradient

        Matrix3 const F_unimodular = unimodular(F);

        auto const nonaffine_stretch = compute_nonaffine_stretch(F_unimodular);

        // Compute the non-affine stretch derivatives
        Matrix3 const h = compute_h_tensor(F_unimodular);
        Matrix6 const H = compute_H_tensor(F_unimodular);

        // Compute the microstress and micro moduli
        auto const micro_kirchhoff_f = G_eff * pade_first(nonaffine_stretch, N) * nonaffine_stretch;

        auto const micro_moduli_f = G_eff * pade_second(nonaffine_stretch, N);

        // Compute the macrostress and macromoduli for chain force
        Matrix3 const macro_kirchhoff_f = micro_kirchhoff_f * std::pow(nonaffine_stretch, 1.0 - p)
                                          * h;

        Matrix6 const macro_moduli_f = (micro_moduli_f * std::pow(nonaffine_stretch, 2.0 - 2.0 * p)
                                        - (p - 1.0) * micro_kirchhoff_f
                                              * std::pow(nonaffine_stretch, 1.0 - 2.0 * p))
                                           * outer_product(h, h)
                                       + micro_kirchhoff_f * std::pow(nonaffine_stretch, 1.0 - p) * H;

        // Compute the non-affine tube contribution
        Matrix3 const k = compute_k_tensor(F_unimodular);

        Matrix6 const K = compute_K_tensor(F_unimodular);
        Matrix6 const G = compute_G_tensor(F_unimodular);

        // Compute the macrostress and macromoduli for tube contraint
        Matrix3 const macro_kirchhoff_c = -G_eff * N * effective_tube_geometry * k;
        Matrix6 const macro_moduli_c = G_eff * N * effective_tube_geometry * (K + G);

        // Superimposed stress response from tube and chain contributions
        Matrix3 const macro_kirchhoff = macro_kirchhoff_f + macro_kirchhoff_c;

        Matrix6 const macro_moduli = macro_moduli_f + macro_moduli_c;

        auto const pressure = J * volumetric_free_energy_dJ(J, K_eff);

        // Perform the deviatoric projection for the stress and macro moduli
        cauchy_stresses[l] = compute_kirchhoff_stress(pressure, macro_kirchhoff) / J;
        tangent_operators[l] = compute_material_tangent(J, K_eff, macro_moduli, macro_kirchhoff);
    }
}

double NonAffineMicrosphere::compute_nonaffine_stretch(Matrix3 const& F_unimodular) const
{
    auto const p = non_affine_stretch_parameter;

    return std::pow(unit_sphere.integrate(0.0,
                                          [&](auto const& coordinates, auto const& l) {
                                              auto const & [ r, r_outer_r ] = coordinates;

                                              Vector3 const t = deformed_tangent(F_unimodular, r);

                                              return std::pow(compute_area_stretch(t), p);
                                          }),
                    1.0 / p);
}

Matrix3 NonAffineMicrosphere::compute_h_tensor(Matrix3 const& F_unimodular) const
{
    auto const p = non_affine_stretch_parameter;

    return unit_sphere.integrate(Matrix3::Zero().eval(), [&](auto const& xyz, auto const& l) -> Matrix3 {
        auto const & [ r, r_o_r ] = xyz;

        Vector3 const t = deformed_tangent(F_unimodular, r);

        return std::pow(compute_microstretch(t), p - 2.0) * outer_product(t, t);
    });
}

Matrix6 NonAffineMicrosphere::compute_H_tensor(Matrix3 const& F_unimodular) const
{
    auto const p = non_affine_stretch_parameter;

    // clang-format off
    return (p - 2.0)
           * unit_sphere.integrate(Matrix6::Zero().eval(), [&](auto const& xyz, auto const& l) -> Matrix6 {
                 auto const & [ r, r_o_r ] = xyz;

                 Vector3 const t = deformed_tangent(F_unimodular, r);

                 return std::pow(compute_microstretch(t), p - 4.0) * outer_product(t, t, t, t);
             });
    // clang-format on
}

Matrix3 NonAffineMicrosphere::compute_k_tensor(Matrix3 const& F_unimodular) const
{
    auto const q = non_affine_tube_parameter;

    // clang-format off
    return q * unit_sphere.integrate(Matrix3::Zero().eval(), [&](auto const& xyz, auto const& l) -> Matrix3 {
        auto const & [ r, r_o_r ] = xyz;

        Vector3 const n = deformed_normal(F_unimodular, r);

        return std::pow(compute_area_stretch(n), q - 2.0) * outer_product(n, n);
    });
    // clang-format on
}

Matrix6 NonAffineMicrosphere::compute_K_tensor(Matrix3 const& F_unimodular) const
{
    auto const q = non_affine_tube_parameter;

    // clang-format off
    return q * (q - 2.0)
           * unit_sphere.integrate(Matrix6::Zero().eval(), [&](auto const& xyz, auto const& l) -> Matrix6 {
                 auto const & [ r, r_o_r ] = xyz;

                 Vector3 const n = deformed_normal(F_unimodular, r);

                 return std::pow(compute_area_stretch(n), q - 4.0) * outer_product(n, n, n, n);
             });
    // clang-format on
}

Matrix6 NonAffineMicrosphere::compute_G_tensor(Matrix3 const& F_unimodular) const
{
    auto const q = non_affine_tube_parameter;

    // clang-format off
    return 2.0 * q
           * unit_sphere.integrate(Matrix6::Zero().eval(), [&](auto const& xyz, auto const& l) -> Matrix6 {
                 auto const & [ r, r_o_r ] = xyz;

                 Vector3 const n = deformed_normal(F_unimodular, r);

                 return std::pow(compute_area_stretch(n), q - 2.0) * compute_o_dot_product(n);
             });
    // clang-format on
}
}
