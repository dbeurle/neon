
#pragma once

#include "constitutive/mechanical/solid/gaussian_affine_microsphere.hpp"

#include "numeric/dense_matrix.hpp"

#include <vector>

namespace neon::mechanical::solid
{
/// \ingroup Hyperelastic
/// \addtogroup Hyperelastic
/// \{
///
/// gaussian_ageing_affine_microsphere is responsible for computing the Cauchy
/// stress and the material tangent in implicit methods with chemical ageing.
///
/// The affine microsphere model is used to model the mechanical response using
/// micromechanical motivations and homogenises the force from a single chain
/// over a unit sphere using a Gaussian statistical description.
///
/// Warning: This method is highly experimental and is a proof-of-concept
/// implementation.
class gaussian_ageing_affine_microsphere : public gaussian_affine_microsphere
{
public:
    /// \param variables Reference to internal state variable store
    /// \param material_data Json object with input file material data
    explicit gaussian_ageing_affine_microsphere(std::shared_ptr<internal_variables_t>& variables,
                                                json const& material_data,
                                                unit_sphere_quadrature::point const rule);

    virtual ~gaussian_ageing_affine_microsphere() = default;

    virtual void update_internal_variables(double const time_step_size) override;

private:
    [[nodiscard]] matrix3 compute_macro_stress(matrix3 const& F_bar,
                                               std::vector<double> const& modulus1,
                                               std::vector<double> const& modulus2,
                                               double const active_segments,
                                               double const reduction_factor) const;

    [[nodiscard]] matrix6 compute_macro_moduli(matrix3 const& F_bar,
                                               std::vector<double> const& modulus2,
                                               std::vector<double> const& modulus3,
                                               double const active_segments,
                                               double const reduction_factor) const;

    /// Compute the force contribution to the secondary modulus
    /// \return evaluated integrand
    [[nodiscard]] double evaluate_integrand_first(double const creation_rate,
                                                  double const reduction_factor,
                                                  double const micro_stretch) const;

    /// Compute the energy contribution to the secondary modulus
    /// \return evaluated integrand
    [[nodiscard]] double evaluate_integrand_second(double const creation_rate,
                                                   double const reduction_factor,
                                                   double const micro_stretch) const;

    /// Compute the energy contribution to the secondary modulus
    /// \return evaluated integrand
    [[nodiscard]] double evaluate_integrand_third(double const creation_rate,
                                                  double const reduction_factor,
                                                  double const micro_stretch) const;

    /// Compute the integral prefactor that involves the free energy and the constant contribution
    /// from the statistical mechanical chain formation
    /// \return factor value
    [[nodiscard]] double compute_prefactor_first(double const micro_stretch,
                                                 double const active_segments) const;

    /// Compute the integral prefactor that involves the free energy and the constant contribution
    /// from the statistical mechanical chain formation
    /// \return factor value
    [[nodiscard]] double compute_prefactor_second(double const micro_stretch,
                                                  double const active_segments) const;

    /// Integrate the time history and return the integrated value
    [[nodiscard]] double integrate_history(double const modulus,
                                           double const sphere_old,
                                           double const sphere,
                                           double const time_step_size) const;

private:
    /// Material with micromechanical parameters
    ageing_micromechanical_elastomer material;

    double last_time_step_size = 0.0;
};
/** \} */
}
