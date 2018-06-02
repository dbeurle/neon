
#include "gaussian_affine_microsphere.hpp"

#include "constitutive/internal_variables.hpp"
#include "constitutive/mechanical/detail/microsphere.hpp"
#include "constitutive/mechanical/volumetric_free_energy.hpp"

#include <tbb/parallel_for.h>

namespace neon::mechanical::solid
{
gaussian_affine_microsphere::gaussian_affine_microsphere(std::shared_ptr<internal_variables_t>& variables,
                                                         json const& material_data,
                                                         unit_sphere_quadrature::Rule const rule)
    : constitutive_model(variables), unit_sphere(rule), material(material_data)
{
    variables->add(internal_variables_t::fourth::tangent_operator);

    // Commit these to history in case of failure on first time step
    variables->commit();
}

void gaussian_affine_microsphere::update_internal_variables(double const time_step_size)
{
    auto& tangent_operators = variables->get(internal_variables_t::fourth::tangent_operator);

    auto const& deformation_gradients = variables->get(
        internal_variables_t::second::DeformationGradient);
    auto& cauchy_stresses = variables->get(internal_variables_t::second::cauchy_stress);

    auto const& det_deformation_gradients = variables->get(internal_variables_t::scalar::DetF);

    auto const K = material.bulk_modulus();
    auto const G = material.shear_modulus();
    auto const N = material.segments_per_chain();

    tbb::parallel_for(std::size_t{0}, deformation_gradients.size(), [&](auto const l) {
        // Compute the deviatoric (unimodular) deformation gradient
        matrix3 const F_bar = unimodular(deformation_gradients[l]);

        auto const J = det_deformation_gradients[l];

        auto const pressure = J * volumetric_free_energy_dJ(J, K);

        // Project the stresses to obtain the Kirchhoff macro-stress
        matrix3 const macro_stress = compute_macro_stress(F_bar, G);

        cauchy_stresses[l] = compute_kirchhoff_stress(pressure, macro_stress) / J;

        tangent_operators[l] = compute_material_tangent(J,
                                                        K,
                                                        compute_macro_moduli(F_bar, G),
                                                        macro_stress);
    });
}

matrix3 gaussian_affine_microsphere::compute_kirchhoff_stress(double const pressure,
                                                              matrix3 const& macro_stress) const
{
    using namespace voigt;

    return pressure * matrix3::Identity() + kinetic::from(P * kinetic::to(macro_stress));
}

matrix6 gaussian_affine_microsphere::compute_material_tangent(double const J,
                                                              double const K,
                                                              matrix6 const& macro_C,
                                                              matrix3 const& macro_stress) const
{
    auto const pressure = J * volumetric_free_energy_dJ(J, K);
    auto const kappa = std::pow(J, 2) * volumetric_free_energy_second_d2J(J, K);

    // clang-format off
    matrix6 const D = macro_C
                    + 2.0 / 3.0 * macro_stress.trace() * voigt::kinematic::identity()
                    - 2.0 / 3.0 * (outer_product(macro_stress, matrix3::Identity()) +
                                   outer_product(matrix3::Identity(), macro_stress));

    // clang-format on
    return (kappa + pressure) * IoI - 2.0 * pressure * I + P * D * P;
}

matrix3 gaussian_affine_microsphere::compute_macro_stress(matrix3 const& F_unimodular,
                                                          double const shear_modulus) const
{
    return 3.0 * shear_modulus
           * unit_sphere.integrate(matrix3::Zero().eval(),
                                   [&](auto const& coordinates, auto const l) -> matrix3 {
                                       auto const& [r, _] = coordinates;

                                       vector3 const t = deformed_tangent(F_unimodular, r);

                                       return outer_product(t, t);
                                   });
}

matrix6 gaussian_affine_microsphere::compute_macro_moduli(matrix3 const& F_unimodular,
                                                          double const shear_modulus) const
{
    return -3.0 * shear_modulus
           * unit_sphere.integrate(matrix6::Zero().eval(),
                                   [&](auto const& coordinates, auto const l) -> matrix6 {
                                       auto const& [r, _] = coordinates;

                                       vector3 const t = deformed_tangent(F_unimodular, r);

                                       auto const micro_stretch = compute_microstretch(t);

                                       return std::pow(micro_stretch, -2) * outer_product(t, t, t, t);
                                   });
}
}
