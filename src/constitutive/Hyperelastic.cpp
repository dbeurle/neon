
#include "Hyperelastic.hpp"

#include "InternalVariables.hpp"
#include "numeric/DenseTypes.hpp"

#include <algorithm>
#include <json/json.h>
#include <range/v3/algorithm.hpp>
#include <range/v3/view.hpp>

namespace neon
{
Hyperelastic::Hyperelastic(InternalVariables& variables) : ConstitutiveModel(variables) {}

NeoHooke::NeoHooke(InternalVariables& variables, Json::Value const& material_data)
    : Hyperelastic(variables), material(material_data)
{
    // The Neo-Hookean model requires the following deformation measures
    // - C, right Cauchy-Green tensor for S, the PK1 stress (total)
    // - B, left Cauchy-Green tensor for kappa, the Kirchoff stress (updated)
    // Which can be computed from the deformation gradient
    variables.add(InternalVariables::Tensor::Kirchhoff,
                  InternalVariables::Tensor::DeformationGradient);

    variables.add(InternalVariables::Matrix::MaterialTangent, 6);
}

void NeoHooke::update_internal_variables()
{
    // Get references into the hash table
    auto[Fvec, Kvec] = variables(InternalVariables::Tensor::DeformationGradient,
                                 InternalVariables::Tensor::Kirchhoff);

    // Perform the update of the internal variables
    std::transform(Fvec.begin(), Fvec.end(), Kvec.begin(), [this](auto const& F) {

        auto const[μ0, λ0] = material.LameConstants();

        // Left Cauchy Green deformation tensor
        Matrix3 B = F.matrix() * F.matrix().transpose();

        auto const J = F.determinant();

        Matrix3 const I = Matrix3::Identity();

        return λ0 * std::log(J) * I + μ0 * (B - I);
    });
}

void NeoHooke::update_continuum_tangent()
{
    auto& DVec = variables(InternalVariables::Matrix::MaterialTangent);
    auto const& Fvec = variables(InternalVariables::Tensor::DeformationGradient);

    using namespace ranges;

    for_each(view::zip(Fvec, DVec), [&](auto const& zip_pair) {

        auto & [ F, D ] = zip_pair;

        // FIXME Move this outside lambda when compiler is fixed
        auto const[μ0, λ0] = material.LameConstants();

        auto const J = F.determinant();

        auto const μ = μ0 - λ0 * std::log(J);

        D << λ0 + 2.0 * μ, λ0, λ0, 0.0, 0.0, 0.0, //
            λ0, λ0 + 2.0 * μ, λ0, 0.0, 0.0, 0.0,  //
            λ0, λ0, λ0 + 2.0 * μ, 0.0, 0.0, 0.0,  //
            0.0, 0.0, 0.0, 2.0 * μ, 0.0, 0.0,     //
            0.0, 0.0, 0.0, 0.0, 2.0 * μ, 0.0,     //
            0.0, 0.0, 0.0, 0.0, 0.0, 2.0 * μ;
    });
}
}
