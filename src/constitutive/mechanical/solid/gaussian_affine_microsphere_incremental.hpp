
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
<<<<<<< HEAD

protected:
    /**
     * Compute the Kirchhoff stress using the deviatoric projection of the
     * macro stress tensor according to
     * \f{align*}{
       \boldsymbol{\tau} &= p \boldsymbol{g}^{-1} + \mathbb{P} : \bar{\boldsymbol{\tau}} \f}
     * @param pressure Hydrostatic pressure
     * @param macro_stress Stress tensor from unit sphere homogenisation
     */
    [[nodiscard]] matrix3 compute_kirchhoff_stress(double const pressure,
                                                   matrix3 const& macro_stress) const;

    /**
     * Compute the material tangent including the sdeviatoric projection defined as
     *\f{align*}{
        \boldsymbol{C} &= (\kappa + p) \mathbf{g}^{-1} \otimes \mathbf{g}^{-1} - 2p\mathbb{I} +
     \mathbb{P} : \left[\bar{\mathbb{C}} + \frac{2}{3}(\boldsymbol{\tau} : \boldsymbol{g})
     \mathbb{I} - \frac{2}{3}
     (\bar{\boldsymbol{\tau}} \otimes \boldsymbol{g}^{-1} + \boldsymbol{g}^{-1} \otimes
     \bar{\boldsymbol{\tau}}) \right] : \mathbb{P} \f}
     * @param jacobian_determinant Determinant of deformation gradient
     * @param shear_modulus Shear modulus
     * @param macro_C Macromoduli from unit sphere
     * @param macro_stress Macrostress from unit sphere
     */
    [[nodiscard]] matrix6 compute_material_tangent(double const jacobian_determinant,
                                                   double const shear_modulus,
                                                   matrix6 const& macro_C,
                                                   matrix3 const& macro_stress) const;

    /**
     * Compute the macro stress using the unit sphere homogenisation
     * technique for a given F and N and perform the deviatoric projection
         \f{align*}{
            \bar{\boldsymbol{\tau}}_f &= 3\mu
                                         \sum_{i=1}^{m}\boldsymbol{\Delta t}_i
                                               \otimes \boldsymbol{\Delta t}_i w_i
         \f}
     * @param F_unimodular Unimodular decomposition of the deformation gradient
     * @param shear_modulus The material shear modulus
     * @param N number of segments per chain
     * @return Kirchhoff stress tensor
     */
    [[nodiscard]] matrix3 compute_macro_stress(matrix3 const& F_unimodular,
                                               double const shear_modulus,
                                               double const N) const;

    /**
     * Compute the material tangent matrix using the unit sphere homogenisation
     * technique for a given F and N
         \f{align*}{
            \bar{\mathbb{C}}_f &= -3\mu \sum_{i=1}^{m} \bar{\lambda}_i^{-2}
                                     \boldsymbol{\Delta t}_i \otimes
                                     \boldsymbol{\Delta t}_i \otimes
                                     \boldsymbol{\Delta t}_i \otimes
                                     \boldsymbol{\Delta t}_i w_i
         \f}
     * @param F_unimodular Unimodular decomposition of the deformation gradient
     * @param bulk_modulus The material bulk modulus
     * @param N number of segments per chain
     * @return Macromoduli from unit sphere homogenisation
     */
    [[nodiscard]] matrix6 compute_macro_moduli(matrix3 const& F_unimodular,
                                               double const bulk_modulus,
                                               double const N) const;
=======
>>>>>>> Remove code inherited from non-incremental version
};
/** \} */
}
