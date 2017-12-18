
#include "IsotropicLinearElasticity.hpp"

#include "Exceptions.hpp"
#include "constitutive/InternalVariables.hpp"

#include "numeric/mechanics"

#include <range/v3/view/transform.hpp>

namespace neon::mechanical::plane
{
IsotropicLinearElasticity::IsotropicLinearElasticity(std::shared_ptr<InternalVariables>& variables,
                                                     Json::Value const& material_data,
                                                     plane const state)
    : ConstitutiveModel(variables), material(material_data), state(state)
{
    variables->add(InternalVariables::Tensor::LinearisedStrain,
                   InternalVariables::Scalar::VonMisesStress);

    // Add material tangent with the linear elasticity spatial moduli
    variables->add(InternalVariables::Matrix::TangentOperator, elastic_moduli());

    std::cout << "Has linearised strain "
              << variables->has(InternalVariables::Tensor::LinearisedStrain) << std::endl;

    std::cout << "Has VonMisesStress strain "
              << variables->has(InternalVariables::Scalar::VonMisesStress) << std::endl;
}

IsotropicLinearElasticity::~IsotropicLinearElasticity() = default;

void IsotropicLinearElasticity::update_internal_variables(double const time_step_size)
{
    using namespace ranges;

    std::cout << "Performing the internal variable updates" << std::endl;

    std::cout << "Has variables " << variables->has(InternalVariables::Tensor::LinearisedStrain)
              << std::endl;

    // Extract the internal variables
    auto [elastic_strains,
          cauchy_stresses] = variables->fetch(InternalVariables::Tensor::LinearisedStrain,
                                              InternalVariables::Tensor::Cauchy);

    std::cout << "Size of elastic strains: " << elastic_strains.size() << std::endl;

    std::cout << "Accessing the von Mises stress" << std::endl;

    auto& von_mises_stresses = variables->fetch(InternalVariables::Scalar::VonMisesStress);

    std::cout << "Size of displacement gradient: "
              << variables->fetch(InternalVariables::Tensor::DisplacementGradient).size()
              << std::endl;

    std::cout << "Computing elastic strains" << std::endl;

    // Compute the linear strain gradient from the displacement gradient
    elastic_strains = variables->fetch(InternalVariables::Tensor::DisplacementGradient)
                      | view::transform([](auto const& H) { return 0.5 * (H + H.transpose()); });

    std::cout << "Computing Cauchy stress" << std::endl;

    // Compute Cauchy stress from the linear elastic strains
    cauchy_stresses = elastic_strains | view::transform([this](auto const& elastic_strain) {
                          return compute_cauchy_stress(elastic_strain);
                      });

    std::cout << "Computing von Mises stress" << std::endl;

    // Compute the von Mises equivalent stress
    von_mises_stresses = cauchy_stresses | view::transform([](auto const& cauchy_stress) {
                             return von_mises_stress(cauchy_stress);
                         });

    std::cout << "All done!" << std::endl;
}

matrix3 IsotropicLinearElasticity::elastic_moduli() const
{
    auto [lambda, shear_modulus] = material.Lame_parameters();

    if (state == plane::stress)
    {
        lambda = 2.0 * lambda * shear_modulus / (lambda + 2.0 * shear_modulus);
    }
    // clang-format off
    return (matrix3() << lambda + 2.0 * shear_modulus, lambda, 0.0,
                         lambda, lambda + 2.0 * shear_modulus, 0.0,
                         0.0, 0.0, shear_modulus).finished();
    // clang-format on
}

matrix2 IsotropicLinearElasticity::compute_cauchy_stress(matrix2 const& elastic_strain) const
{
    auto const G = material.shear_modulus();
    auto const lambda_e = material.lambda();
    return lambda_e * elastic_strain.trace() * matrix2::Identity() + 2.0 * G * elastic_strain;
}
}
