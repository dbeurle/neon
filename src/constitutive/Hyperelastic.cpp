
#include "Hyperelastic.hpp"

#include "InternalVariables.hpp"
#include "numeric/DenseTypes.hpp"

#include <algorithm>
#include <json/json.h>
#include <range/v3/algorithm.hpp>
#include <range/v3/view.hpp>

#include <iostream>

namespace neon
{
Hyperelastic::Hyperelastic(InternalVariables& variables) : ConstitutiveModel(variables) {}

NeoHooke::NeoHooke(InternalVariables& variables, Json::Value const& material_data)
    : Hyperelastic(variables), material(material_data)
{
    // The Neo-Hookean model requires the following deformation measures
    // - C, right Cauchy-Green tensor for S, the PK1 stress (total)
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

    σ = view::zip(F, detF) | view::transform([this, &I](auto const& FdetF) -> Matrix3 {
            auto const[μ0, λ0] = material.LameConstants();

            auto const & [ F, J ] = FdetF;

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
}
