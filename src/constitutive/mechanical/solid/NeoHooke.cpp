
#include "NeoHooke.hpp"

#include "constitutive/InternalVariables.hpp"

#include <range/v3/algorithm/for_each.hpp>
#include <range/v3/view/transform.hpp>
#include <range/v3/view/zip.hpp>

namespace neon::mechanical::solid
{
NeoHooke::NeoHooke(std::shared_ptr<InternalVariables>& variables, Json::Value const& material_data)
    : ConstitutiveModel(variables), material(material_data)
{
    // The Neo-Hookean model requires the deformation gradient and the Cauchy
    // stress, which are both allocated by default in the mesh object
    variables->add(InternalVariables::rank4::tangent_operator);
}

void NeoHooke::update_internal_variables(double const time_step_size)
{
    using namespace ranges;

    // Get references into the hash table
    auto [F_list, cauchy_stresses] = variables->fetch(InternalVariables::Tensor::DeformationGradient,
                                                      InternalVariables::Tensor::Cauchy);

    auto& tangent_operators = variables->fetch(InternalVariables::rank4::tangent_operator);
    auto const& detF_list = variables->fetch(InternalVariables::Scalar::DetF);

    auto const I = matrix3::Identity();

    // Compute stresses
    cauchy_stresses = view::zip(F_list, detF_list)
                      | view::transform([this, &I](auto const& tpl) -> matrix3 {
                            auto const [lambda, shear_modulus] = material.Lame_parameters();

                            auto const& [F, J] = tpl;

                            // Left Cauchy Green deformation tensor
                            matrix3 const B = F * F.transpose();

                            // Compute Kirchhoff stress and transform to Cauchy
                            return (lambda * std::log(J) * I + shear_modulus * (B - I)) / J;
                        });

    // Compute tangent moduli
    for_each(view::zip(tangent_operators, detF_list), [&](auto const& tpl) {
        auto& [D, J] = tpl;

        auto const [lambda, shear_modulus_0] = material.Lame_parameters();

        auto const shear_modulus = shear_modulus_0 - lambda * std::log(J);

        D << lambda + 2.0 * shear_modulus, lambda, lambda, 0.0, 0.0, 0.0, //
            lambda, lambda + 2.0 * shear_modulus, lambda, 0.0, 0.0, 0.0,  //
            lambda, lambda, lambda + 2.0 * shear_modulus, 0.0, 0.0, 0.0,  //
            0.0, 0.0, 0.0, shear_modulus, 0.0, 0.0,                       //
            0.0, 0.0, 0.0, 0.0, shear_modulus, 0.0,                       //
            0.0, 0.0, 0.0, 0.0, 0.0, shear_modulus;
    });
}
}
