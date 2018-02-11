
#include "isotropic_linear_elasticity.hpp"

#include "Exceptions.hpp"
#include "constitutive/internal_variables.hpp"

#include "numeric/mechanics"

#include <range/v3/view/transform.hpp>

namespace neon::mechanical::solid
{
isotropic_linear_elasticity::isotropic_linear_elasticity(std::shared_ptr<internal_variables_t>& variables,
                                                     json const& material_data)
    : constitutive_model(variables), material(material_data)
{
    variables->add(internal_variables_t::Tensor::LinearisedStrain,
                   internal_variables_t::Scalar::VonMisesStress);

    // Add material tangent with the linear elasticity spatial moduli
    variables->add(internal_variables_t::rank4::tangent_operator, elastic_moduli());
}

isotropic_linear_elasticity::~isotropic_linear_elasticity() = default;

void isotropic_linear_elasticity::update_internal_variables(double const time_step_size)
{
    using namespace ranges;

    // Extract the internal variables
    auto [elastic_strains,
          cauchy_stresses] = variables->fetch(internal_variables_t::Tensor::LinearisedStrain,
                                              internal_variables_t::Tensor::Cauchy);

    auto& von_mises_stresses = variables->fetch(internal_variables_t::Scalar::VonMisesStress);

    // Compute the linear strain gradient from the displacement gradient
    elastic_strains = variables->fetch(internal_variables_t::Tensor::DisplacementGradient)
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

matrix6 isotropic_linear_elasticity::elastic_moduli() const
{
    auto const [lambda, shear_modulus] = material.Lame_parameters();

    // clang-format off
    return (matrix6() << lambda + 2.0 * shear_modulus, lambda, lambda, 0.0, 0.0, 0.0,
                         lambda, lambda + 2.0 * shear_modulus, lambda, 0.0, 0.0, 0.0,
                         lambda, lambda, lambda + 2.0 * shear_modulus, 0.0, 0.0, 0.0,
                         0.0, 0.0, 0.0, shear_modulus, 0.0, 0.0,
                         0.0, 0.0, 0.0, 0.0, shear_modulus, 0.0,
                         0.0, 0.0, 0.0, 0.0, 0.0, shear_modulus).finished();
    // clang-format on
}
}
