
#include "IsotropicLinearElasticity.hpp"

#include "Exceptions.hpp"
#include "constitutive/InternalVariables.hpp"

#include <range/v3/view/transform.hpp>

namespace neon::mech::plane
{
IsotropicLinearElasticity::IsotropicLinearElasticity(InternalVariables& variables,
                                                     Json::Value const& material_data,
                                                     State const state)
    : ConstitutiveModel(variables), material(material_data), state(state)
{
    variables.add(InternalVariables::Tensor::LinearisedStrain);
    variables.add(InternalVariables::Scalar::VonMisesStress);

    // Add material tangent with the linear elasticity spatial moduli
    variables.add(InternalVariables::Matrix::TangentOperator, elastic_moduli());
}

IsotropicLinearElasticity::~IsotropicLinearElasticity() = default;

void IsotropicLinearElasticity::update_internal_variables(double const time_step_size)
{
    using namespace ranges;

    // Extract the internal variables
    auto[elastic_strains, cauchy_stresses] = variables(InternalVariables::Tensor::LinearisedStrain,
                                                       InternalVariables::Tensor::Cauchy);

    auto& von_mises_stresses = variables(InternalVariables::Scalar::VonMisesStress);

    // Compute the linear strain gradient from the displacement gradient
    elastic_strains = variables(InternalVariables::Tensor::DisplacementGradient)
                      | view::transform([](auto const& H) { return 0.5 * (H + H.transpose()); });

    // Compute Cauchy stress from the linear elastic strains
    cauchy_stresses = elastic_strains | view::transform([this](auto const& elastic_strain) {
                          return compute_cauchy_stress(elastic_strain);
                      });

    // Compute the von Mises equivalent stress
    von_mises_stresses = cauchy_stresses | view::transform([](auto const& cauchy_stress) {
                             return von_mises_stress(cauchy_stress);
                         });
}

Matrix3 IsotropicLinearElasticity::elastic_moduli() const
{
    auto[lambda, shear_modulus] = material.Lame_parameters();

    if (state == State::PlaneStress)
    {
        lambda = 2.0 * lambda * shear_modulus / (lambda + 2.0 * shear_modulus);
    }
    // clang-format off
    return (Matrix3() << lambda + 2.0 * shear_modulus, lambda, 0.0,
                         lambda, lambda + 2.0 * shear_modulus, 0.0,
                         0.0, 0.0, shear_modulus).finished();
    // clang-format on
}

Matrix2 IsotropicLinearElasticity::compute_cauchy_stress(Matrix2 const& elastic_strain) const
{
    auto const G = material.shear_modulus();
    auto const lambda_e = material.lambda();
    return lambda_e * elastic_strain.trace() * Matrix2::Identity() + 2.0 * G * elastic_strain;
}
}
