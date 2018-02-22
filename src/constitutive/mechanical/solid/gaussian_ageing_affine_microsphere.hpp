
#pragma once

#include "constitutive/mechanical/solid/gaussian_affine_microsphere.hpp"

/**
 * \ingroup Hyperelastic
 * \addtogroup Hyperelastic
 * \{
 *
 * gaussian_ageing_affine_microsphere is responsible for computing the Cauchy
 * stress and the material tangent in implicit methods when ageing is present.
 *
 * The affine microsphere model is used to model elastomer materials using
 * micromechanical motivations and homogenises the force from a single chain
 * over a unit sphere using a Gaussian statistical description.
 *
 * Warning: This method is highly experimental and this is a proof-of-concept
 * implementation.
 */
class gaussian_ageing_affine_microsphere : public gaussian_affine_microsphere
{
public:
    /**
     * @param variables Reference to internal state variable store
     * @param material_data Json object with input file material data
     */
    explicit gaussian_ageing_affine_microsphere(std::shared_ptr<internal_variables_t>& variables,
                                                json const& material_data,
                                                unit_sphere_quadrature::Rule const rule);

    virtual void update_internal_variables(double const time_step_size) override;

private:
    stochastic_micromechanical_elastomer material; //!< Material with micromechanical parameters
};
/** \} */
