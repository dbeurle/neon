
#include "IsotropicLinearElasticity.hpp"

#include "Exceptions.hpp"
#include "constitutive/InternalVariables.hpp"

#include "numeric/mechanics"

#include <range/v3/view/transform.hpp>
#include <range/v3/view/zip.hpp>

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
}

IsotropicLinearElasticity::~IsotropicLinearElasticity() = default;

void IsotropicLinearElasticity::update_internal_variables(double const time_step_size)
{
    using namespace ranges;

    // std::cout << "Performing the internal variable updates" << std::endl;

    // Extract the internal variables
    auto [elastic_strains,
          cauchy_stresses] = variables->fetch(InternalVariables::Tensor::LinearisedStrain,
                                              InternalVariables::Tensor::Cauchy);

    auto& von_mises_stresses = variables->fetch(InternalVariables::Scalar::VonMisesStress);

    auto const& tangents = variables->fetch(InternalVariables::Matrix::TangentOperator);

    // std::cout << "Computing elastic strains" << std::endl;

    // Compute the linear strain gradient from the displacement gradient
    elastic_strains = variables->fetch(InternalVariables::Tensor::DisplacementGradient)
                      | view::transform([](auto const& H) { return 0.5 * (H + H.transpose()); });

    for (auto const& elastic_strain : elastic_strains)
    {
        if (elastic_strain.hasNaN()) std::cout << "NAN DETECTED\n";
    }

    // Compute Cauchy stress from the linear elastic strains
    cauchy_stresses = view::zip(tangents, elastic_strains) | view::transform([this](auto const& tpl) {
                          auto const& [C, elastic_strain] = tpl;
                          return voigt::kinetic::from(C * voigt::kinematic::to(elastic_strain));
                      });

    for (auto const& cauchy_stress : cauchy_stresses)
    {
        if (cauchy_stress.hasNaN()) std::cout << "NAN DETECTED\n";
    }

    // for (auto const& cauchy_stress : cauchy_stresses)
    // {
    //     std::cout << cauchy_stress << std::endl << std::endl;
    // }

    // Compute the von Mises equivalent stress
    von_mises_stresses = cauchy_stresses | view::transform([](auto const& cauchy_stress) {
                             return von_mises_stress(cauchy_stress);
                         });

    // std::cout << "All done!" << std::endl;
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
