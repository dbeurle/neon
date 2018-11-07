
#include "compressible_neohooke.hpp"

#include "constitutive/internal_variables.hpp"

namespace neon::mechanics::solid
{
compressible_neohooke::compressible_neohooke(std::shared_ptr<internal_variables_t>& variables,
                                             json const& material_data)
    : constitutive_model(variables), material(material_data)
{
    // The Neo-Hookean model requires the deformation gradient and the Cauchy
    // stress, which are both allocated by default in the mesh object
    variables->add(variable::fourth::tangent_operator);
}

void compressible_neohooke::update_internal_variables(double const time_step_size)
{
    auto [deformation_gradients,
          cauchy_stresses] = variables->get(variable::second::deformation_gradient,
                                            variable::second::cauchy_stress);

    auto& tangent_operators = variables->get(variable::fourth::tangent_operator);
    auto const& determinants = variables->get(variable::scalar::DetF);

    // compute stresses
    std::transform(begin(deformation_gradients),
                   end(deformation_gradients),
                   begin(determinants),
                   begin(cauchy_stresses),
                   [&](matrix3 const& F, double const J) -> matrix3 {
                       auto const [lambda, shear_modulus] = material.Lame_parameters();

                       auto const I = matrix3::Identity();

                       // Left Cauchy Green deformation tensor
                       matrix3 const B = F * F.transpose();

                       // Compute Kirchhoff stress and transform to Cauchy
                       return (lambda * std::log(J) * I + shear_modulus * (B - I)) / J;
                   });

    // compute material tangent operators
    std::transform(begin(determinants),
                   end(determinants),
                   begin(tangent_operators),
                   [&](double const J) -> matrix6 {
                       auto const [lambda, shear_modulus_0] = material.Lame_parameters();

                       auto const shear_modulus = shear_modulus_0 - lambda * std::log(J);

                       return (matrix6() << lambda + 2.0 * shear_modulus,
                               lambda,
                               lambda,
                               0.0,
                               0.0,
                               0.0, //
                               lambda,
                               lambda + 2.0 * shear_modulus,
                               lambda,
                               0.0,
                               0.0,
                               0.0, //
                               lambda,
                               lambda,
                               lambda + 2.0 * shear_modulus,
                               0.0,
                               0.0,
                               0.0, //
                               0.0,
                               0.0,
                               0.0,
                               shear_modulus,
                               0.0,
                               0.0, //
                               0.0,
                               0.0,
                               0.0,
                               0.0,
                               shear_modulus,
                               0.0, //
                               0.0,
                               0.0,
                               0.0,
                               0.0,
                               0.0,
                               shear_modulus)
                           .finished();
                   });
}
}
