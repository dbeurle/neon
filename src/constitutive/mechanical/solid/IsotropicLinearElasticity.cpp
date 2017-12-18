
#include "IsotropicLinearElasticity.hpp"

#include "Exceptions.hpp"
#include "constitutive/InternalVariables.hpp"

#include "numeric/mechanics"

#include <range/v3/view/transform.hpp>

namespace neon::mechanical::solid
{
IsotropicLinearElasticity::IsotropicLinearElasticity(std::shared_ptr<InternalVariables>& variables,
                                                     Json::Value const& material_data)
    : ConstitutiveModel(variables), material(material_data)
{
    variables->add(InternalVariables::Tensor::LinearisedStrain,
                   InternalVariables::Scalar::VonMisesStress);

    // Add material tangent with the linear elasticity spatial moduli
    variables->add(InternalVariables::Matrix::TangentOperator, elastic_moduli());
}

IsotropicLinearElasticity::~IsotropicLinearElasticity() = default;

void IsotropicLinearElasticity::update_internal_variables(double const time_step_size)
{
    using namespace ranges;

    // Extract the internal variables
    auto [elastic_strains,
          cauchy_stresses] = variables->fetch(InternalVariables::Tensor::LinearisedStrain,
                                              InternalVariables::Tensor::Cauchy);

    auto& von_mises_stresses = variables->fetch(InternalVariables::Scalar::VonMisesStress);

    // Compute the linear strain gradient from the displacement gradient
    elastic_strains = variables->fetch(InternalVariables::Tensor::DisplacementGradient)
                      | view::transform([](auto const& H) { return 0.5 * (H + H.transpose()); });

    // Compute Cauchy stress from the linear elastic strains
    cauchy_stresses = elastic_strains | view::transform([this](auto const& elastic_strain) {
                          return compute_cauchy_stress(material.shear_modulus(),
                                                       material.lambda(),
                                                       elastic_strain);
                      });

    // Compute the von Mises equivalent stress
    von_mises_stresses = cauchy_stresses | view::transform([](auto const& cauchy_stress) {
                             return von_mises_stress(cauchy_stress);
                         });
}

Matrix6 IsotropicLinearElasticity::elastic_moduli() const
{
    auto const [lambda, shear_modulus] = material.Lame_parameters();

    // clang-format off
    return (Matrix6() << lambda + 2.0 * shear_modulus, lambda, lambda, 0.0, 0.0, 0.0,
                         lambda, lambda + 2.0 * shear_modulus, lambda, 0.0, 0.0, 0.0,
                         lambda, lambda, lambda + 2.0 * shear_modulus, 0.0, 0.0, 0.0,
                         0.0, 0.0, 0.0, shear_modulus, 0.0, 0.0,
                         0.0, 0.0, 0.0, 0.0, shear_modulus, 0.0,
                         0.0, 0.0, 0.0, 0.0, 0.0, shear_modulus).finished();
    // clang-format on
}
}
