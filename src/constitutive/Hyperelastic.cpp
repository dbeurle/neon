
#include "Hyperelastic.hpp"

#include "InternalVariables.hpp"
#include "numeric/DenseTypes.hpp"

#include <algorithm>
#include <json/json.h>
#include <range/v3/algorithm.hpp>
#include <range/v3/view.hpp>

#include <exception>
#include <iostream>

namespace neon
{
Hyperelastic::Hyperelastic(InternalVariables& variables) : ConstitutiveModel(variables) {}

NeoHooke::NeoHooke(InternalVariables& variables, Json::Value const& material_data)
    : Hyperelastic(variables), material(material_data)
{
    // The Neo-Hookean model requires the following deformation measures
    // - B, left Cauchy-Green tensor for kappa, the Kirchoff stress (updated)
    // Which can be computed from the deformation gradient and transformed
    // to the Cauchy stress with very little effort
    variables.add(InternalVariables::Matrix::MaterialTangent, 6);
}

void NeoHooke::update_internal_variables()
{
    using namespace ranges;

    // Get references into the hash table
    auto[F, σ] = variables(InternalVariables::Tensor::DeformationGradient,
                           InternalVariables::Tensor::CauchyStress);

    auto const& detF = variables(InternalVariables::Scalar::DetF);

    auto const I = Matrix3::Identity();

    σ = view::zip(F, detF) | view::transform([this, &I](auto const& tpl) -> Matrix3 {
            auto const[μ0, λ0] = material.LameConstants();

            auto const & [ F, J ] = tpl;

            // Left Cauchy Green deformation tensor
            auto const B = F * F.transpose();

            return (λ0 * std::log(J) * I + μ0 * (B - I)) / J;
        });
}

void NeoHooke::update_continuum_tangent()
{
    auto& D_list = variables(InternalVariables::Matrix::MaterialTangent);
    auto const& F_list = variables(InternalVariables::Tensor::DeformationGradient);
    auto const& detF_list = variables(InternalVariables::Scalar::DetF);

    ranges::for_each(ranges::view::zip(F_list, D_list, detF_list), [&](auto const& FDdetF) {

        auto & [ F, D, J ] = FDdetF;

        auto const[μ0, λ0] = material.LameConstants();

        auto const μ = μ0 - λ0 * std::log(J);

        D << λ0 + 2.0 * μ, λ0, λ0, 0.0, 0.0, 0.0, //
            λ0, λ0 + 2.0 * μ, λ0, 0.0, 0.0, 0.0,  //
            λ0, λ0, λ0 + 2.0 * μ, 0.0, 0.0, 0.0,  //
            0.0, 0.0, 0.0, μ, 0.0, 0.0,           //
            0.0, 0.0, 0.0, 0.0, μ, 0.0,           //
            0.0, 0.0, 0.0, 0.0, 0.0, μ;
    });
}

AffineMicrosphere::AffineMicrosphere(InternalVariables& variables, Json::Value const& material_data)
    : Hyperelastic(variables), material(material_data)
{
    if (material_data["SegmentsPerChain"].empty())
        throw std::runtime_error("SegmentsPerChain not specified in material data\n");

    segments_per_chain = material_data["SegmentsPerChain"].asInt();

    variables.add(InternalVariables::Matrix::MaterialTangent, 6);

    // Deviatoric stress
    variables.add(InternalVariables::Tensor::Kirchhoff);
}

void AffineMicrosphere::update_internal_variables()
{
    using namespace ranges;

    // Get references into the hash table
    auto[F_list, σ_list, τ_list] = variables(InternalVariables::Tensor::DeformationGradient,
                                             InternalVariables::Tensor::CauchyStress,
                                             InternalVariables::Tensor::Kirchhoff);

    auto const& detF_list = variables(InternalVariables::Scalar::DetF);

    auto N = segments_per_chain;

    τ_list = view::zip(F_list, detF_list) | view::transform([this, &N](auto const& tpl) -> Matrix3 {
                 auto const[μ0, λ0] = material.LameConstants();

                 auto const & [ F, J ] = tpl;

                 // Unimodular decomposition of F
                 Matrix3 const unimodular_F = std::pow(J, -1.0 / 3.0) * F;

                 return μ0 * unit_sphere.integrate(Matrix3::Zero(),
                                                   [&N, &unimodular_F](auto const& coordinates,
                                                                       auto const& l) -> Matrix3 {
                                                       auto const & [ r, r_outer_r ] = coordinates;

                                                       // Deformed tangents
                                                       auto const t = unimodular_F * r;

                                                       // Microstretches
                                                       auto const λ = t.norm();

                                                       return (3.0 * N - std::pow(λ, 2)) /
                                                              (N - std::pow(λ, 2)) * t *
                                                              t.transpose();
                                                   });
             });

    // Perform the projection of the stresses
    σ_list = view::zip(τ_list, detF_list) | view::transform([this](auto const& τdetF) -> Matrix3 {

                 auto const[μ0, λ0] = material.LameConstants();

                 auto const & [ τ_dev, J ] = τdetF;

                 auto const pressure = J * volumetric_free_energy_derivative(J, μ0);

                 return deviatoric_projection(pressure, τ_dev) / J;
             });
}

void AffineMicrosphere::update_continuum_tangent()
{
    using namespace ranges;

    auto& D_list = variables(InternalVariables::Matrix::MaterialTangent);

    auto const & [ F_list, τ_list ] = variables(InternalVariables::Tensor::DeformationGradient,
                                                InternalVariables::Tensor::Kirchhoff);

    auto const& detF_list = variables(InternalVariables::Scalar::DetF);

    auto N = segments_per_chain;

    for_each(view::zip(F_list, τ_list, D_list, detF_list), [&N, this](auto const& tpl) {

        auto const[μ0, λ0] = material.LameConstants();

        auto & [ F, τ_dev, D, J ] = tpl;

        auto const pressure = J * volumetric_free_energy_derivative(J, μ0);
        auto const κ = std::pow(J, 2) * volumetric_free_energy_second_derivative(J, μ0);

        // Unimodular decomposition of F
        Matrix3 const unimodular_F = std::pow(J, -1.0 / 3.0) * F;

        Matrix const IoI = I_outer_I();
        Matrix const I = fourth_order_identity();

        auto const D_dev =
            deviatoric_projection(unit_sphere
                                      .integrate(Matrix::Zero(6, 6),
                                                 [&](auto const& coordinates, auto const& l) -> Matrix {
                                                     auto const & [ r, r_outer_r ] = coordinates;

                                                     // Deformed tangents
                                                     auto const t = unimodular_F * r;

                                                     // Microstretches
                                                     auto const λ = t.norm();

                                                     auto const a =
                                                         std::pow(λ, -2) *
                                                         ((std::pow(λ, 4) + 3.0 * std::pow(N, 2)) /
                                                              std::pow(N - std::pow(λ, 2), 2) -
                                                          (3.0 * N - std::pow(λ, 2)) /
                                                              (N - std::pow(λ, 2)));

                                                     return a * t_outer_t_outer_t_outer_t(t);
                                                 }),
                                  τ_dev);

        auto const D_vol = (κ + pressure) * IoI - 2.0 * pressure * I;

        D = D_dev + D_vol;
    });
}

double AffineMicrosphere::volumetric_free_energy_derivative(double const J,
                                                            double const bulk_modulus) const
{
    return bulk_modulus / 2.0 * (J - 1.0 / J);
}
double AffineMicrosphere::volumetric_free_energy_second_derivative(double const J,
                                                                   double const bulk_modulus) const
{
    return bulk_modulus / 2.0 * (1.0 + 1.0 / std::pow(J, 2));
}

Matrix AffineMicrosphere::t_outer_t_outer_t_outer_t(Vector3 const& t) const
{
    Matrix v(6, 6);
    v << std::pow(t(0), 4),                    //
        std::pow(t(0), 2) * std::pow(t(1), 2), //
        std::pow(t(0), 2) * std::pow(t(2), 2), //
        std::pow(t(0), 2) * t(1) * t(2),       //
        std::pow(t(0), 3) * t(2),              //
        std::pow(t(0), 3) * t(1),              //

        std::pow(t(0), 2) * std::pow(t(1), 2), //
        std::pow(t(1), 4),                     //
        std::pow(t(1), 2) * std::pow(t(2), 2), //
        std::pow(t(1), 3) * t(2),              //
        t(0) * std::pow(t(1), 2) * t(2),       //
        t(0) * std::pow(t(1), 3),              //

        std::pow(t(0), 2) * std::pow(t(2), 2), //
        std::pow(t(1), 2) * std::pow(t(2), 2), //
        std::pow(t(2), 4),                     //
        t(1) * std::pow(t(2), 3),              //
        t(0) * std::pow(t(2), 3),              //
        t(0) * t(1) * std::pow(t(2), 2),       //

        std::pow(t(0), 2) * t(1) * t(2),       //
        std::pow(t(1), 3) * t(2),              //
        t(1) * std::pow(t(2), 3),              //
        std::pow(t(1), 2) * std::pow(t(2), 2), //
        t(0) * t(1) * std::pow(t(2), 2),       //
        t(0) * std::pow(t(1), 2) * t(2),       //

        std::pow(t(0), 3) * t(2),              //
        t(0) * std::pow(t(1), 2) * t(2),       //
        t(0) * std::pow(t(2), 3),              //
        t(0) * t(1) * std::pow(t(2), 2),       //
        std::pow(t(0), 2) * std::pow(t(2), 2), //
        std::pow(t(0), 2) * t(1) * t(2),       //

        std::pow(t(0), 3) * t(1),        //
        t(0) * std::pow(t(1), 3),        //
        t(0) * t(1) * std::pow(t(2), 2), //
        t(0) * std::pow(t(1), 2) * t(2), //
        std::pow(t(0), 2) * t(1) * t(2), //
        std::pow(t(0), 2) * std::pow(t(1), 2);
    return v;
}

Matrix3 AffineMicrosphere::deviatoric_projection(double const pressure, Matrix3 const& τ_dev) const
{
    Matrix3 P_double_dot_τ;
    P_double_dot_τ << 2 * τ_dev(0, 0) / 3.0 - τ_dev(1, 1) / 3.0 - τ_dev(2, 2) / 3.0, //
        τ_dev(0, 1) / 2.0 + τ_dev(1, 0) / 2.0,                                       //
        τ_dev(0, 2) / 2.0 + τ_dev(2, 0) / 2.0,                                       //
                                                                                     //
        τ_dev(0, 1) / 2.0 + τ_dev(1, 0) / 2.0,                                       //
        -τ_dev(0, 0) / 3.0 + 2 * τ_dev(1, 1) / 3.0 - τ_dev(2, 2) / 3.0,              //
        τ_dev(1, 2) / 2.0 + τ_dev(2, 1) / 2.0,                                       //
                                                                                     //
        τ_dev(0, 2) / 2.0 + τ_dev(2, 0) / 2.0,                                       //
        τ_dev(1, 2) / 2.0 + τ_dev(2, 1) / 2.0,                                       //
        -τ_dev(0, 0) / 3.0 - τ_dev(1, 1) / 3.0 + 2 * τ_dev(2, 2) / 3.0;
    return pressure * Matrix3::Identity() + P_double_dot_τ;
}

Matrix AffineMicrosphere::deviatoric_projection(Matrix const& C_dev, Matrix3 const& τ_dev) const
{
    Matrix C(6, 6);
    C << -2 * (C_dev(0, 1) - 2 * (τ_dev(0, 0) + τ_dev(1, 1)) / 3.0) / 9.0 -
             2 * (C_dev(0, 2) - 2 * (τ_dev(0, 0) + τ_dev(2, 2)) / 3.0) / 9.0 -
             2 * (C_dev(1, 0) - 2 * (τ_dev(0, 0) + τ_dev(1, 1)) / 3.0) / 9.0 +
             (C_dev(1, 2) - 2 * (τ_dev(1, 1) + τ_dev(2, 2)) / 3.0) / 9.0 -
             2 * (C_dev(2, 0) - 2 * (τ_dev(0, 0) + τ_dev(2, 2)) / 3.0) / 9.0 +
             (C_dev(2, 1) - 2 * (τ_dev(1, 1) + τ_dev(2, 2)) / 3.0) / 9.0 +
             4 * (C_dev(0, 0) + 2 * τ_dev.trace() / 3.0 - 4 * τ_dev(0, 0) / 3.0) / 9.0 +
             (C_dev(1, 1) + 2 * τ_dev.trace() / 3.0 - 4 * τ_dev(1, 1) / 3.0) / 9.0 +
             (C_dev(2, 2) + 2 * τ_dev.trace() / 3.0 - 4 * τ_dev(2, 2) / 3.0) / 9.0, //
        4 * (C_dev(0, 1) - 2 * (τ_dev(0, 0) + τ_dev(1, 1)) / 3.0) / 9.0 -
            2 * (C_dev(0, 2) - 2 * (τ_dev(0, 0) + τ_dev(2, 2)) / 3.0) / 9.0 +
            (C_dev(1, 0) - 2 * (τ_dev(0, 0) + τ_dev(1, 1)) / 3.0) / 9.0 +
            (C_dev(1, 2) - 2 * (τ_dev(1, 1) + τ_dev(2, 2)) / 3.0) / 9.0 +
            (C_dev(2, 0) - 2 * (τ_dev(0, 0) + τ_dev(2, 2)) / 3.0) / 9.0 -
            2 * (C_dev(2, 1) - 2 * (τ_dev(1, 1) + τ_dev(2, 2)) / 3.0) / 9.0 -
            2 * (C_dev(0, 0) + 2 * τ_dev.trace() / 3.0 - 4 * τ_dev(0, 0) / 3.0) / 9.0 -
            2 * (C_dev(1, 1) + 2 * τ_dev.trace() / 3.0 - 4 * τ_dev(1, 1) / 3.0) / 9.0 +
            (C_dev(2, 2) + 2 * τ_dev.trace() / 3.0 - 4 * τ_dev(2, 2) / 3.0) / 9.0, //
        -2 * (C_dev(0, 1) - 2 * (τ_dev(0, 0) + τ_dev(1, 1)) / 3.0) / 9.0 +
            4 * (C_dev(0, 2) - 2 * (τ_dev(0, 0) + τ_dev(2, 2)) / 3.0) / 9.0 +
            (C_dev(1, 0) - 2 * (τ_dev(0, 0) + τ_dev(1, 1)) / 3.0) / 9.0 -
            2 * (C_dev(1, 2) - 2 * (τ_dev(1, 1) + τ_dev(2, 2)) / 3.0) / 9.0 +
            (C_dev(2, 0) - 2 * (τ_dev(0, 0) + τ_dev(2, 2)) / 3.0) / 9.0 +
            (C_dev(2, 1) - 2 * (τ_dev(1, 1) + τ_dev(2, 2)) / 3.0) / 9.0 -
            2 * (C_dev(0, 0) + 2 * τ_dev.trace() / 3.0 - 4 * τ_dev(0, 0) / 3.0) / 9.0 +
            (C_dev(1, 1) + 2 * τ_dev.trace() / 3.0 - 4 * τ_dev(1, 1) / 3.0) / 9.0 -
            2 * (C_dev(2, 2) + 2 * τ_dev.trace() / 3.0 - 4 * τ_dev(2, 2) / 3.0) / 9.0, //
        (C_dev(0, 3) - 2 * τ_dev(1, 2) / 3.0) / 3.0 + (C_dev(0, 3) - 2 * τ_dev(2, 1) / 3.0) / 3.0 -
            (C_dev(1, 3) - 2 * τ_dev(1, 2) / 3.0) / 6.0 -
            (C_dev(1, 3) - 2 * τ_dev(2, 1) / 3.0) / 6.0 -
            (C_dev(2, 3) - 2 * τ_dev(1, 2) / 3.0) / 6.0 -
            (C_dev(2, 3) - 2 * τ_dev(2, 1) / 3.0) / 6.0, //
        (C_dev(0, 4) - 2 * τ_dev(0, 2) / 3.0) / 3.0 + (C_dev(0, 4) - 2 * τ_dev(2, 0) / 3.0) / 3.0 -
            (C_dev(1, 4) - 2 * τ_dev(0, 2) / 3.0) / 6.0 -
            (C_dev(1, 4) - 2 * τ_dev(2, 0) / 3.0) / 6.0 -
            (C_dev(2, 4) - 2 * τ_dev(0, 2) / 3.0) / 6.0 - (C_dev(2, 4) - 2 * τ_dev(2, 0) / 3.0) / 6.0,
        (C_dev(0, 5) - 2 * τ_dev(0, 1) / 3.0) / 3.0 + (C_dev(0, 5) - 2 * τ_dev(1, 0) / 3.0) / 3.0 -
            (C_dev(1, 5) - 2 * τ_dev(0, 1) / 3.0) / 6.0 -
            (C_dev(1, 5) - 2 * τ_dev(1, 0) / 3.0) / 6.0 -
            (C_dev(2, 5) - 2 * τ_dev(0, 1) / 3.0) / 6.0 -
            (C_dev(2, 5) - 2 * τ_dev(1, 0) / 3.0) / 6.0, //
        //
        (C_dev(0, 1) - 2 * (τ_dev(0, 0) + τ_dev(1, 1)) / 3.0) / 9.0 +
            (C_dev(0, 2) - 2 * (τ_dev(0, 0) + τ_dev(2, 2)) / 3.0) / 9.0 +
            4 * (C_dev(1, 0) - 2 * (τ_dev(0, 0) + τ_dev(1, 1)) / 3.0) / 9.0 -
            2 * (C_dev(1, 2) - 2 * (τ_dev(1, 1) + τ_dev(2, 2)) / 3.0) / 9.0 -
            2 * (C_dev(2, 0) - 2 * (τ_dev(0, 0) + τ_dev(2, 2)) / 3.0) / 9.0 +
            (C_dev(2, 1) - 2 * (τ_dev(1, 1) + τ_dev(2, 2)) / 3.0) / 9.0 -
            2 * (C_dev(0, 0) + 2 * τ_dev.trace() / 3.0 - 4 * τ_dev(0, 0) / 3.0) / 9.0 -
            2 * (C_dev(1, 1) + 2 * τ_dev.trace() / 3.0 - 4 * τ_dev(1, 1) / 3.0) / 9.0 +
            (C_dev(2, 2) + 2 * τ_dev.trace() / 3.0 - 4 * τ_dev(2, 2) / 3.0) / 9.0, //
        -2 * (C_dev(0, 1) - 2 * (τ_dev(0, 0) + τ_dev(1, 1)) / 3.0) / 9.0 +
            (C_dev(0, 2) - 2 * (τ_dev(0, 0) + τ_dev(2, 2)) / 3.0) / 9.0 -
            2 * (C_dev(1, 0) - 2 * (τ_dev(0, 0) + τ_dev(1, 1)) / 3.0) / 9.0 -
            2 * (C_dev(1, 2) - 2 * (τ_dev(1, 1) + τ_dev(2, 2)) / 3.0) / 9.0 +
            (C_dev(2, 0) - 2 * (τ_dev(0, 0) + τ_dev(2, 2)) / 3.0) / 9.0 -
            2 * (C_dev(2, 1) - 2 * (τ_dev(1, 1) + τ_dev(2, 2)) / 3.0) / 9.0 +
            (C_dev(0, 0) + 2 * τ_dev.trace() / 3.0 - 4 * τ_dev(0, 0) / 3.0) / 9.0 +
            4 * (C_dev(1, 1) + 2 * τ_dev.trace() / 3.0 - 4 * τ_dev(1, 1) / 3.0) / 9.0 +
            (C_dev(2, 2) + 2 * τ_dev.trace() / 3.0 - 4 * τ_dev(2, 2) / 3.0) / 9.0, //
        (C_dev(0, 1) - 2 * (τ_dev(0, 0) + τ_dev(1, 1)) / 3.0) / 9.0 -
            2 * (C_dev(0, 2) - 2 * (τ_dev(0, 0) + τ_dev(2, 2)) / 3.0) / 9.0 -
            2 * (C_dev(1, 0) - 2 * (τ_dev(0, 0) + τ_dev(1, 1)) / 3.0) / 9.0 +
            4 * (C_dev(1, 2) - 2 * (τ_dev(1, 1) + τ_dev(2, 2)) / 3.0) / 9.0 +
            (C_dev(2, 0) - 2 * (τ_dev(0, 0) + τ_dev(2, 2)) / 3.0) / 9.0 +
            (C_dev(2, 1) - 2 * (τ_dev(1, 1) + τ_dev(2, 2)) / 3.0) / 9.0 +
            (C_dev(0, 0) + 2 * τ_dev.trace() / 3.0 - 4 * τ_dev(0, 0) / 3.0) / 9.0 -
            2 * (C_dev(1, 1) + 2 * τ_dev.trace() / 3.0 - 4 * τ_dev(1, 1) / 3.0) / 9.0 -
            2 * (C_dev(2, 2) + 2 * τ_dev.trace() / 3.0 - 4 * τ_dev(2, 2) / 3.0) / 9.0,
        -(C_dev(0, 3) - 2 * τ_dev(1, 2) / 3.0) / 6.0 - (C_dev(0, 3) - 2 * τ_dev(2, 1) / 3.0) / 6.0 +
            (C_dev(1, 3) - 2 * τ_dev(1, 2) / 3.0) / 3.0 +
            (C_dev(1, 3) - 2 * τ_dev(2, 1) / 3.0) / 3.0 -
            (C_dev(2, 3) - 2 * τ_dev(1, 2) / 3.0) / 6.0 -
            (C_dev(2, 3) - 2 * τ_dev(2, 1) / 3.0) / 6.0, //
        -(C_dev(0, 4) - 2 * τ_dev(0, 2) / 3.0) / 6.0 - (C_dev(0, 4) - 2 * τ_dev(2, 0) / 3.0) / 6.0 +
            (C_dev(1, 4) - 2 * τ_dev(0, 2) / 3.0) / 3.0 +
            (C_dev(1, 4) - 2 * τ_dev(2, 0) / 3.0) / 3.0 -
            (C_dev(2, 4) - 2 * τ_dev(0, 2) / 3.0) / 6.0 -
            (C_dev(2, 4) - 2 * τ_dev(2, 0) / 3.0) / 6.0, //
        -(C_dev(0, 5) - 2 * τ_dev(0, 1) / 3.0) / 6.0 - (C_dev(0, 5) - 2 * τ_dev(1, 0) / 3.0) / 6.0 +
            (C_dev(1, 5) - 2 * τ_dev(0, 1) / 3.0) / 3.0 +
            (C_dev(1, 5) - 2 * τ_dev(1, 0) / 3.0) / 3.0 -
            (C_dev(2, 5) - 2 * τ_dev(0, 1) / 3.0) / 6.0 -
            (C_dev(2, 5) - 2 * τ_dev(1, 0) / 3.0) / 6.0, //
        //
        (C_dev(0, 1) - 2 * (τ_dev(0, 0) + τ_dev(1, 1)) / 3.0) / 9.0 +
            (C_dev(0, 2) - 2 * (τ_dev(0, 0) + τ_dev(2, 2)) / 3.0) / 9.0 -
            2 * (C_dev(1, 0) - 2 * (τ_dev(0, 0) + τ_dev(1, 1)) / 3.0) / 9.0 +
            (C_dev(1, 2) - 2 * (τ_dev(1, 1) + τ_dev(2, 2)) / 3.0) / 9.0 +
            4 * (C_dev(2, 0) - 2 * (τ_dev(0, 0) + τ_dev(2, 2)) / 3.0) / 9.0 -
            2 * (C_dev(2, 1) - 2 * (τ_dev(1, 1) + τ_dev(2, 2)) / 3.0) / 9.0 -
            2 * (C_dev(0, 0) + 2 * τ_dev.trace() / 3.0 - 4 * τ_dev(0, 0) / 3.0) / 9.0 +
            (C_dev(1, 1) + 2 * τ_dev.trace() / 3.0 - 4 * τ_dev(1, 1) / 3.0) / 9.0 -
            2 * (C_dev(2, 2) + 2 * τ_dev.trace() / 3.0 - 4 * τ_dev(2, 2) / 3.0) / 9.0, //
        -2 * (C_dev(0, 1) - 2 * (τ_dev(0, 0) + τ_dev(1, 1)) / 3.0) / 9.0 +
            (C_dev(0, 2) - 2 * (τ_dev(0, 0) + τ_dev(2, 2)) / 3.0) / 9.0 +
            (C_dev(1, 0) - 2 * (τ_dev(0, 0) + τ_dev(1, 1)) / 3.0) / 9.0 +
            (C_dev(1, 2) - 2 * (τ_dev(1, 1) + τ_dev(2, 2)) / 3.0) / 9.0 -
            2 * (C_dev(2, 0) - 2 * (τ_dev(0, 0) + τ_dev(2, 2)) / 3.0) / 9.0 +
            4 * (C_dev(2, 1) - 2 * (τ_dev(1, 1) + τ_dev(2, 2)) / 3.0) / 9.0 +
            (C_dev(0, 0) + 2 * τ_dev.trace() / 3.0 - 4 * τ_dev(0, 0) / 3.0) / 9.0 -
            2 * (C_dev(1, 1) + 2 * τ_dev.trace() / 3.0 - 4 * τ_dev(1, 1) / 3.0) / 9.0 -
            2 * (C_dev(2, 2) + 2 * τ_dev.trace() / 3.0 - 4 * τ_dev(2, 2) / 3.0) / 9.0, //
        (C_dev(0, 1) - 2 * (τ_dev(0, 0) + τ_dev(1, 1)) / 3.0) / 9.0 -
            2 * (C_dev(0, 2) - 2 * (τ_dev(0, 0) + τ_dev(2, 2)) / 3.0) / 9.0 +
            (C_dev(1, 0) - 2 * (τ_dev(0, 0) + τ_dev(1, 1)) / 3.0) / 9.0 -
            2 * (C_dev(1, 2) - 2 * (τ_dev(1, 1) + τ_dev(2, 2)) / 3.0) / 9.0 -
            2 * (C_dev(2, 0) - 2 * (τ_dev(0, 0) + τ_dev(2, 2)) / 3.0) / 9.0 -
            2 * (C_dev(2, 1) - 2 * (τ_dev(1, 1) + τ_dev(2, 2)) / 3.0) / 9.0 +
            (C_dev(0, 0) + 2 * τ_dev.trace() / 3.0 - 4 * τ_dev(0, 0) / 3.0) / 9.0 +
            (C_dev(1, 1) + 2 * τ_dev.trace() / 3.0 - 4 * τ_dev(1, 1) / 3.0) / 9.0 +
            4 * (C_dev(2, 2) + 2 * τ_dev.trace() / 3.0 - 4 * τ_dev(2, 2) / 3.0) / 9.0,
        -(C_dev(0, 3) - 2 * τ_dev(1, 2) / 3.0) / 6.0 - (C_dev(0, 3) - 2 * τ_dev(2, 1) / 3.0) / 6.0 -
            (C_dev(1, 3) - 2 * τ_dev(1, 2) / 3.0) / 6.0 -
            (C_dev(1, 3) - 2 * τ_dev(2, 1) / 3.0) / 6.0 +
            (C_dev(2, 3) - 2 * τ_dev(1, 2) / 3.0) / 3.0 +
            (C_dev(2, 3) - 2 * τ_dev(2, 1) / 3.0) / 3.0, //
        -(C_dev(0, 4) - 2 * τ_dev(0, 2) / 3.0) / 6.0 - (C_dev(0, 4) - 2 * τ_dev(2, 0) / 3.0) / 6.0 -
            (C_dev(1, 4) - 2 * τ_dev(0, 2) / 3.0) / 6.0 -
            (C_dev(1, 4) - 2 * τ_dev(2, 0) / 3.0) / 6.0 +
            (C_dev(2, 4) - 2 * τ_dev(0, 2) / 3.0) / 3.0 + (C_dev(2, 4) - 2 * τ_dev(2, 0) / 3.0) / 3.0,
        -(C_dev(0, 5) - 2 * τ_dev(0, 1) / 3.0) / 6.0 - (C_dev(0, 5) - 2 * τ_dev(1, 0) / 3.0) / 6.0 -
            (C_dev(1, 5) - 2 * τ_dev(0, 1) / 3.0) / 6.0 -
            (C_dev(1, 5) - 2 * τ_dev(1, 0) / 3.0) / 6.0 +
            (C_dev(2, 5) - 2 * τ_dev(0, 1) / 3.0) / 3.0 +
            (C_dev(2, 5) - 2 * τ_dev(1, 0) / 3.0) / 3.0, //
        //
        (C_dev(3, 0) - 2 * τ_dev(1, 2) / 3.0) / 3.0 + (C_dev(3, 0) - 2 * τ_dev(2, 1) / 3.0) / 3.0 -
            (C_dev(3, 1) - 2 * τ_dev(1, 2) / 3.0) / 6.0 -
            (C_dev(3, 1) - 2 * τ_dev(2, 1) / 3.0) / 6.0 -
            (C_dev(3, 2) - 2 * τ_dev(1, 2) / 3.0) / 6.0 - (C_dev(3, 2) - 2 * τ_dev(2, 1) / 3.0) / 6.0,
        -(C_dev(3, 0) - 2 * τ_dev(1, 2) / 3.0) / 6.0 - (C_dev(3, 0) - 2 * τ_dev(2, 1) / 3.0) / 6.0 +
            (C_dev(3, 1) - 2 * τ_dev(1, 2) / 3.0) / 3.0 +
            (C_dev(3, 1) - 2 * τ_dev(2, 1) / 3.0) / 3.0 -
            (C_dev(3, 2) - 2 * τ_dev(1, 2) / 3.0) / 6.0 - (C_dev(3, 2) - 2 * τ_dev(2, 1) / 3.0) / 6.0,
        -(C_dev(3, 0) - 2 * τ_dev(1, 2) / 3.0) / 6.0 - (C_dev(3, 0) - 2 * τ_dev(2, 1) / 3.0) / 6.0 -
            (C_dev(3, 1) - 2 * τ_dev(1, 2) / 3.0) / 6.0 -
            (C_dev(3, 1) - 2 * τ_dev(2, 1) / 3.0) / 6.0 +
            (C_dev(3, 2) - 2 * τ_dev(1, 2) / 3.0) / 3.0 +
            (C_dev(3, 2) - 2 * τ_dev(2, 1) / 3.0) / 3.0,             //
        C_dev(3, 3) + τ_dev.trace() / 3.0, C_dev(3, 4), C_dev(3, 5), //
        //
        (C_dev(4, 0) - 2 * τ_dev(0, 2) / 3.0) / 3.0 + (C_dev(4, 0) - 2 * τ_dev(2, 0) / 3.0) / 3.0 -
            (C_dev(4, 1) - 2 * τ_dev(0, 2) / 3.0) / 6.0 -
            (C_dev(4, 1) - 2 * τ_dev(2, 0) / 3.0) / 6.0 -
            (C_dev(4, 2) - 2 * τ_dev(0, 2) / 3.0) / 6.0 -
            (C_dev(4, 2) - 2 * τ_dev(2, 0) / 3.0) / 6.0, //
        -(C_dev(4, 0) - 2 * τ_dev(0, 2) / 3.0) / 6.0 - (C_dev(4, 0) - 2 * τ_dev(2, 0) / 3.0) / 6.0 +
            (C_dev(4, 1) - 2 * τ_dev(0, 2) / 3.0) / 3.0 +
            (C_dev(4, 1) - 2 * τ_dev(2, 0) / 3.0) / 3.0 -
            (C_dev(4, 2) - 2 * τ_dev(0, 2) / 3.0) / 6.0 -
            (C_dev(4, 2) - 2 * τ_dev(2, 0) / 3.0) / 6.0, //
        -(C_dev(4, 0) - 2 * τ_dev(0, 2) / 3.0) / 6.0 - (C_dev(4, 0) - 2 * τ_dev(2, 0) / 3.0) / 6.0 -
            (C_dev(4, 1) - 2 * τ_dev(0, 2) / 3.0) / 6.0 -
            (C_dev(4, 1) - 2 * τ_dev(2, 0) / 3.0) / 6.0 +
            (C_dev(4, 2) - 2 * τ_dev(0, 2) / 3.0) / 3.0 +
            (C_dev(4, 2) - 2 * τ_dev(2, 0) / 3.0) / 3.0,             //
        C_dev(4, 3), C_dev(4, 4) + τ_dev.trace() / 3.0, C_dev(4, 5), //
        //
        (C_dev(5, 0) - 2 * τ_dev(0, 1) / 3.0) / 3.0 + (C_dev(5, 0) - 2 * τ_dev(1, 0) / 3.0) / 3.0 -
            (C_dev(5, 1) - 2 * τ_dev(0, 1) / 3.0) / 6.0 -
            (C_dev(5, 1) - 2 * τ_dev(1, 0) / 3.0) / 6.0 -
            (C_dev(5, 2) - 2 * τ_dev(0, 1) / 3.0) / 6.0 -
            (C_dev(5, 2) - 2 * τ_dev(1, 0) / 3.0) / 6.0, //
        -(C_dev(5, 0) - 2 * τ_dev(0, 1) / 3.0) / 6.0 - (C_dev(5, 0) - 2 * τ_dev(1, 0) / 3.0) / 6.0 +
            (C_dev(5, 1) - 2 * τ_dev(0, 1) / 3.0) / 3.0 +
            (C_dev(5, 1) - 2 * τ_dev(1, 0) / 3.0) / 3.0 -
            (C_dev(5, 2) - 2 * τ_dev(0, 1) / 3.0) / 6.0 -
            (C_dev(5, 2) - 2 * τ_dev(1, 0) / 3.0) / 6.0, //
        -(C_dev(5, 0) - 2 * τ_dev(0, 1) / 3.0) / 6.0 - (C_dev(5, 0) - 2 * τ_dev(1, 0) / 3.0) / 6.0 -
            (C_dev(5, 1) - 2 * τ_dev(0, 1) / 3.0) / 6.0 -
            (C_dev(5, 1) - 2 * τ_dev(1, 0) / 3.0) / 6.0 +
            (C_dev(5, 2) - 2 * τ_dev(0, 1) / 3.0) / 3.0 +
            (C_dev(5, 2) - 2 * τ_dev(1, 0) / 3.0) / 3.0, //
        C_dev(5, 3), C_dev(5, 4), C_dev(5, 5) + τ_dev.trace() / 3.0;

    return C;
}
}
