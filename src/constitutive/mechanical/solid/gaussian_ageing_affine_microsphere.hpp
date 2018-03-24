
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
    /// Store the active and inactive value
    using value_pair = std::pair<double, double>;

    /// Type for the composition of the active/inactive segment history
    using segment_composition = std::vector<value_pair>;

    /// Type for the composition of the active/inactive shear modulus history
    using shear_modulus_composition = std::vector<value_pair>;

    /// Type for the deformation history
    using deformation_composition = std::vector<matrix3>;

public:
    /// \param variables Reference to internal state variable store
    /// \param material_data Json object with input file material data
    explicit gaussian_ageing_affine_microsphere(std::shared_ptr<internal_variables_t>& variables,
                                                json const& material_data,
                                                unit_sphere_quadrature::Rule const rule);

    virtual void update_internal_variables(double const time_step_size) override;

private:
    ageing_micromechanical_elastomer material; /// Material with micromechanical parameters

    std::vector<shear_modulus_composition> shear_moduli; /// Quadrature point shear moduli

    std::vector<segment_composition> segments; /// Quadrature point average segment

    std::vector<deformation_composition> intermediate_deformations; /// Secondary network variables
};
/** \} */
}
