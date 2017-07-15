
#include "NeoHooke.hpp"

#include "InternalVariables.hpp"

#include <range/v3/algorithm.hpp>
#include <range/v3/view.hpp>

namespace neon
{
NeoHooke::NeoHooke(InternalVariables& variables, Json::Value const& material_data)
    : Hyperelastic(variables), material(material_data)
{
    // The Neo-Hookean model requires the following deformation measures
    // - B, left Cauchy-Green tensor for kappa, the Kirchoff stress (updated)
    // Which can be computed from the deformation gradient and transformed
    // to the Cauchy stress with very little effort
    variables.add(InternalVariables::Matrix::TruesdellModuli, 6);
}

void NeoHooke::update_internal_variables(double const Δt)
{
    using namespace ranges;

    // Get references into the hash table
    auto[F_list, σ] =
        variables(InternalVariables::Tensor::DeformationGradient, InternalVariables::Tensor::Cauchy);

    auto& D_list = variables(InternalVariables::Matrix::TruesdellModuli);
    auto const& detF_list = variables(InternalVariables::Scalar::DetF);

    auto const I = Matrix3::Identity();

    // Compute stresses
    σ = view::zip(F_list, detF_list) | view::transform([this, &I](auto const& tpl) -> Matrix3 {
            auto const[λ0, μ0] = material.Lame_parameters();

            auto const & [ F, J ] = tpl;

            // Left Cauchy Green deformation tensor
            auto const B = F * F.transpose();

            return (λ0 * std::log(J) * I + μ0 * (B - I)) / J;
        });

    // Compute tangent moduli
    ranges::for_each(ranges::view::zip(F_list, D_list, detF_list), [&](auto const& tpl) {

        auto & [ F, D, J ] = tpl;

        auto const[λ0, μ0] = material.Lame_parameters();

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
