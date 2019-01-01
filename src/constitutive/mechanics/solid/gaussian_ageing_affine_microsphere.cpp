
#include "constitutive/mechanics/solid/gaussian_ageing_affine_microsphere.hpp"

#include "constitutive/internal_variables.hpp"
#include "constitutive/mechanics/detail/microsphere.hpp"
#include "constitutive/mechanics/volumetric_free_energy.hpp"
#include "quadrature/trapezoidal.hpp"
#include "solver/time/runge_kutta_integration.hpp"
#include "numeric/float_compare.hpp"

#include <tbb/parallel_for.h>

#include <cassert>
#include <numeric>

namespace neon::mechanics::solid
{
gaussian_ageing_affine_microsphere::gaussian_ageing_affine_microsphere(
    std::shared_ptr<internal_variables_t>& variables,
    json const& material_data,
    unit_sphere_quadrature::point const rule)
    : gaussian_affine_microsphere{variables, material_data, rule}, material{material_data}
{
    variables->add(variable::scalar::active_shear_modulus,
                   variable::scalar::inactive_shear_modulus,
                   variable::scalar::active_segments,
                   variable::scalar::inactive_segments,
                   variable::scalar::reduction_factor);

    variables->add(variable::second::intermediate_secondary_kirchhoff_stress);

    names.emplace("active_shear_modulus");
    names.emplace("inactive_shear_modulus");
    names.emplace("active_segments");
    names.emplace("inactive_segments");
    names.emplace("reduction_factor");

    auto [active_shear_modulus,
          active_segments,
          reduction] = variables->get(variable::scalar::active_shear_modulus,
                                      variable::scalar::active_segments,
                                      variable::scalar::reduction_factor);

    std::fill(begin(active_shear_modulus), end(active_shear_modulus), material.shear_modulus());
    std::fill(begin(active_segments), end(active_segments), material.segments_per_chain());
    std::fill(begin(reduction), end(reduction), 1.0);

    variables->commit();
}

void gaussian_ageing_affine_microsphere::update_internal_variables(double const time_step_size)
{
    // Polymer network values
    auto& active_shear_modulus = variables->get(variable::scalar::active_shear_modulus);
    auto& inactive_shear_modulus = variables->get(variable::scalar::inactive_shear_modulus);
    auto& reduction_factor = variables->get(variable::scalar::reduction_factor);
    auto& active_segments = variables->get(variable::scalar::active_segments);
    auto& inactive_segments = variables->get(variable::scalar::inactive_segments);

    auto const& deformation_gradients = variables->get(variable::second::deformation_gradient);
    auto const& deformation_gradients_old = variables->get_old(variable::second::deformation_gradient);

    auto& intermediate_kirchhoff_stresses = variables->get(
        variable::second::intermediate_secondary_kirchhoff_stress);

    auto const& intermediate_kirchhoff_stresses_old = variables->get_old(
        variable::second::intermediate_secondary_kirchhoff_stress);

    auto const& det_F = variables->get(variable::scalar::DetF);

    auto& cauchy_stresses = variables->get(variable::second::cauchy_stress);
    auto& tangent_operators = variables->get(variable::fourth::tangent_operator);

    auto const K{material.bulk_modulus()};

    tbb::parallel_for(std::size_t{0}, deformation_gradients.size(), [&, this](auto const index) {
        // unimodular deformation gradient
        matrix3 const F_bar = unimodular(deformation_gradients[index]);

        // incremental deformation gradient

        if (!is_approx(time_step_size, 0.0))
        {
            // Integrate the network evolution differential equations through
            // the micromechanical material by passing our network parameters
            vector5 parameters;
            parameters(0) = active_shear_modulus[index];
            parameters(1) = inactive_shear_modulus[index];
            parameters(2) = reduction_factor[index];
            parameters(3) = active_segments[index];
            parameters(4) = inactive_segments[index];

            parameters = material.integrate(parameters, time_step_size);

            // Update the history variables for plotting
            active_shear_modulus[index] = parameters(0);
            inactive_shear_modulus[index] = parameters(1);
            reduction_factor[index] = parameters(2);
            active_segments[index] = parameters(3);
            inactive_segments[index] = parameters(4);

            last_time_step_size = time_step_size;
        }

        auto const creation_rate = material.creation_rate(active_shear_modulus[index],
                                                          inactive_shear_modulus[index],
                                                          active_segments[index],
                                                          inactive_segments[index]);

        // explicit shear modulus creation computation (to fix)
        auto const shear_modulus_created = creation_rate * last_time_step_size;

        matrix3 macro_stress = compute_initial_macro_stress(F_bar, reduction_factor[index]);

        matrix3 intermediate_kirchhoff_stress = compute_intermediate_macro_stress(shear_modulus_created,
                                                                                  reduction_factor[index]);

        matrix3 const F_increment = incremental_deformation_gradient(deformation_gradients[index],
                                                                     deformation_gradients_old[index]);

        // push forward old intermediate stress and accumulate with the recently
        // computed intermediate stress
        matrix3 const
            secondary_macro_stress = push_forward_contravariant(F_increment,
                                                                intermediate_kirchhoff_stresses_old[index])
                                     + intermediate_kirchhoff_stress;

        macro_stress += deviatoric(3.0 * reduction_factor[index] * secondary_macro_stress);

        matrix6 const macro_moduli = 0.0
                                     * compute_macro_moduli(F_bar, creation_rate, last_time_step_size);

        auto const J = det_F[index];

        auto const pressure = J * volumetric_free_energy_dJ(J, K);

        intermediate_kirchhoff_stresses[index] = secondary_macro_stress;

        cauchy_stresses[index] = compute_kirchhoff_stress(pressure, macro_stress) / J;

        tangent_operators[index] = compute_material_tangent(J, K, macro_moduli, macro_stress);
    });
}

matrix3 gaussian_ageing_affine_microsphere::compute_initial_macro_stress(matrix3 const& F_bar,
                                                                         double const reduction_factor) const
{
    return 3.0 * reduction_factor * material.shear_modulus()
           * unit_sphere.integrate(matrix3::Zero().eval(),
                                   [&](auto const& coordinates, auto) -> matrix3 {
                                       auto const& [r, _] = coordinates;

                                       vector3 const t = deformed_tangent(F_bar, r);

                                       return outer_product(t, t);
                                   });
}

matrix3 gaussian_ageing_affine_microsphere::compute_intermediate_macro_stress(
    double const shear_modulus_created,
    double const reduction_factor) const
{
    return shear_modulus_created / reduction_factor
           * unit_sphere.integrate(matrix3::Zero().eval(),
                                   [&](auto const& coordinates, auto) -> matrix3 {
                                       auto const& [r, r_outer_r] = coordinates;

                                       return r_outer_r;
                                   });
}

matrix6 gaussian_ageing_affine_microsphere::compute_macro_moduli(matrix3 const& F_bar,
                                                                 double const creation_rate,
                                                                 double const time_step_size) const
{
    return -3.0 / 2.0 * creation_rate * time_step_size
           * unit_sphere.integrate(matrix6::Zero().eval(),
                                   [&](auto const& coordinates, auto) -> matrix6 {
                                       auto const& [r, _] = coordinates;

                                       vector3 const t = deformed_tangent(F_bar, r);

                                       auto const stretch = compute_microstretch(t);

                                       return std::pow(stretch, -3) * outer_product(t, t, t, t);
                                   });
}

}
