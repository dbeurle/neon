
#include "isotropic_linear_elasticity.hpp"

#include "Exceptions.hpp"
#include "constitutive/internal_variables.hpp"

#include "numeric/mechanics"

#include <range/v3/view/transform.hpp>
#include <range/v3/view/zip.hpp>

namespace neon::mechanical::plane
{
isotropic_linear_elasticity::isotropic_linear_elasticity(
    std::shared_ptr<internal_variables_t>& variables, json const& material_data, plane const state)
    : constitutive_model(variables), material(material_data), state(state)
{
    variables->add(internal_variables_t::Tensor::LinearisedStrain,
                   internal_variables_t::scalar::VonMisesStress);

    // Add material tangent with the linear elasticity spatial moduli
    variables->add(internal_variables_t::rank4::tangent_operator, elastic_moduli());
}

isotropic_linear_elasticity::~isotropic_linear_elasticity() = default;

void isotropic_linear_elasticity::update_internal_variables(double const time_step_size)
{
    using namespace ranges;

    // std::cout << "Performing the internal variable updates" << std::endl;

    // Extract the internal variables
    auto [elastic_strains,
          cauchy_stresses] = variables->fetch(internal_variables_t::Tensor::LinearisedStrain,
                                              internal_variables_t::Tensor::CauchyStress);

    auto& von_mises_stresses = variables->fetch(internal_variables_t::scalar::VonMisesStress);

    auto const& tangents = variables->fetch(internal_variables_t::rank4::tangent_operator);

    // std::cout << "Computing elastic strains" << std::endl;

    // Compute the linear strain gradient from the displacement gradient
    elastic_strains = variables->fetch(internal_variables_t::Tensor::DisplacementGradient)
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

    // Compute the von Mises equivalent stress
    von_mises_stresses = cauchy_stresses | view::transform([](auto const& cauchy_stress) {
                             return von_mises_stress(cauchy_stress);
                         });
}

matrix3 isotropic_linear_elasticity::elastic_moduli() const
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
}
