
#pragma once

#include "constitutive/mechanical/solid/gaussian_affine_microsphere.hpp"

namespace neon::mechanical::solid
{
/**
 * \ingroup Hyperelastic
 * \addtogroup Hyperelastic
 * \{
 *
 * gaussian_affine_microsphere_incremental is responsible for computing the
 * Cauchy stress and the material tangent in implicit methods.  The affine
 * microsphere model \cite Miehe2004 is used to model elastomer materials using
 * micromechanical motivations and homogenises the force from a single chain over
 * a unit sphere.  In this variant the non-Gaussian chain theory is replaced
 * with the Gaussian chain theory and made into an incremental step.
 *
 * This constitutive model requires the use of a quadrature scheme for the unit
 * sphere and this internal variable update can be computationally expensive.
 * \sa unit_sphere_quadrature
 */
class gaussian_affine_microsphere_incremental : public gaussian_affine_microsphere
{
public:
    using gaussian_affine_microsphere::gaussian_affine_microsphere;

    virtual void update_internal_variables(double const time_step_size) override;

private:
    MicromechanicalElastomer material; //!< Material with micromechanical parameters
};
/** \} */
}
