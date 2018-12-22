
#pragma once

/// @file

#include "constitutive/constitutive_model.hpp"

#include "material/micromechanical_elastomer.hpp"
#include "numeric/tensor_operations.hpp"
#include "quadrature/sphere/unit_sphere_quadrature.hpp"

#include "io/json_forward.hpp"

namespace neon::mechanics::solid
{
/// \ingroup Hyperelastic
/// \{
/// gaussian_affine_microsphere is responsible for computing the Cauchy stress and the
/// material tangent in implicit methods.  The affine microsphere model @cite Miehe2004
/// is used to model elastomer materials using micromechanical motivations and
/// homogenises the force from a single chain over a unit sphere.  In this variant
/// the non-Gaussian chain theory is replaced with the Gaussian chain theory.
/// This constitutive model requires the use of a quadrature scheme for the unit
/// sphere and this internal variable update can be computationally expensive.
class gaussian_affine_microsphere : public constitutive_model
{
public:
    /// \param variables Reference to internal state variable store
    /// \param material_data json object with input file material data
    explicit gaussian_affine_microsphere(std::shared_ptr<internal_variables_t>& variables,
                                         json const& material_data,
                                         unit_sphere_quadrature::point const p);

    virtual ~gaussian_affine_microsphere() = default;

    virtual void update_internal_variables(double) override;

    material_property const& intrinsic_material() const override final { return material; }

    bool is_finite_deformation() const override final { return true; }

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
     * Compute the material tangent including the deviatoric projection defined as
     *\f{align*}{
        \boldsymbol{C} &= (\kappa + p) \mathbf{g}^{-1} \otimes \mathbf{g}^{-1} - 2p\mathbb{I} +
     \mathbb{P} : \left[\bar{\mathbb{C}} + \frac{2}{3}(\boldsymbol{\tau} : \boldsymbol{g})
     \mathbb{I} - \frac{2}{3}
     (\bar{\boldsymbol{\tau}} \otimes \boldsymbol{g}^{-1} + \boldsymbol{g}^{-1} \otimes
     \bar{\boldsymbol{\tau}}) \right] : \mathbb{P} \f}
     * @param J Determinant of deformation gradient
     * @param K Bulk modulus
     * @param macro_C Macromoduli from unit sphere
     * @param macro_stress Macrostress from unit sphere
     */
    [[nodiscard]] matrix6 compute_material_tangent(double const J,
                                                   double const K,
                                                   matrix6 const& macro_C,
                                                   matrix3 const& macro_stress) const;

    /**
     * Compute the macro stress using the unit sphere homogenisation
     * technique for a given F and N and perform the deviatoric projection
         \f{align*}{
            \bar{\boldsymbol{\tau}}_f &= 3\mu \sum_{i=1}^{m}\boldsymbol{t}_i
                                         \otimes \boldsymbol{t}_i w_i
         \f}
     * @param F_unimodular Unimodular decomposition of the deformation gradient
     * @param shear_modulus The material shear modulus
     * @return Kirchhoff stress tensor
     */
    [[nodiscard]] matrix3 compute_macro_stress(matrix3 const& F_unimodular,
                                               double const shear_modulus) const;

protected:
    /// Unit sphere quadrature rule
    unit_sphere_quadrature unit_sphere;
    /// Outer product
    matrix6 const IoI = voigt::I_outer_I();
    /// Fourth order identity
    matrix6 const I = voigt::kinematic::fourth_order_identity();
    /// Deviatoric fourth order tensor
    matrix6 const P = voigt::kinetic::deviatoric();

    /// Material with micromechanical parameters
    micromechanical_elastomer material;
};
/** \} */
}
