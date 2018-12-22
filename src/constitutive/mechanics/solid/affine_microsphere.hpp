
#pragma once

/// @file

#include "constitutive/constitutive_model.hpp"

#include "material/micromechanical_elastomer.hpp"
#include "numeric/tensor_operations.hpp"
#include "quadrature/sphere/unit_sphere_quadrature.hpp"
#include "io/json_forward.hpp"

namespace neon::mechanics::solid
{
/**
 * \ingroup Hyperelastic
 * \{
 *
 * affine_microsphere is responsible for computing the Cauchy stress and the
 * material tangent in implicit methods.  The affine microsphere model @cite Miehe2004
 * is used to model elastomer materials using micromechanical motivations and
 * homogenises the force from a single chain over a unit sphere.
 *
 * This constitutive model requires the use of a quadrature scheme for the unit
 * sphere and this internal variable update can be computationally expensive and
 * is therefore multithreaded.
 */
class affine_microsphere : public constitutive_model
{
public:
    /// \param variables Reference to internal state variable store
    /// \param material_data Json object with input file material data
    explicit affine_microsphere(std::shared_ptr<internal_variables_t>& variables,
                                json const& material_data,
                                unit_sphere_quadrature::scheme const p);

    virtual ~affine_microsphere() = default;

    void update_internal_variables(double) override;

    [[nodiscard]] material_property const& intrinsic_material() const noexcept override final
    {
        return material;
    }

    [[nodiscard]] bool is_finite_deformation() const noexcept override final { return true; }

protected:
    /**
     * Compute the Kirchhoff stress using the deviatoric projection of the
     * macro stress tensor according to
     * \f{align*}{
       \boldsymbol{\tau} &= p \boldsymbol{g}^{-1} + \mathbb{P} : \bar{\boldsymbol{\tau}} \f}
     * \param pressure Hydrostatic pressure
     * \param macro_stress Stress tensor from unit sphere homogenisation
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
     * \param J Determinant of deformation gradient
     * \param K Shear modulus
     * \param macro_C Macromoduli from unit sphere
     * \param macro_stress Macrostress from unit sphere
     */
    [[nodiscard]] matrix6 compute_material_tangent(double const J,
                                                   double const K,
                                                   matrix6 const& macro_C,
                                                   matrix3 const& macro_stress) const;

    /// Compute the macro stress using the unit sphere homogenisation
    /// technique for a given F and N and perform the deviatoric projection
    /// \param F_unimodular Unimodular decomposition of the deformation gradient
    /// \param shear_modulus The material shear modulus
    /// \param N number of segments per chain
    /// \return Kirchhoff stress tensor
    [[nodiscard]] matrix3 compute_macro_stress(matrix3 const& F_unimodular,
                                               double const shear_modulus,
                                               double const N) const;

    /// Compute the material tangent matrix using the unit sphere homogenisation
    /// technique for a given F and N
    /// \param F_unimodular Unimodular decomposition of the deformation gradient
    /// \param shear_modulus The material shear modulus
    /// \param N number of segments per chain
    /// \return Macromoduli from unit sphere homogenisation
    [[nodiscard]] matrix6 compute_macro_moduli(matrix3 const& F_unimodular,
                                               double const shear_modulus,
                                               double const N) const;

protected:
    /// Unit sphere quadrature rule
    unit_sphere_quadrature unit_sphere;
    /// Outer product
    matrix6 const IoI = voigt::I_outer_I();
    /// Fourth order identity
    matrix6 const I = voigt::kinematic::fourth_order_identity();
    /// Deviatoric fourth order tensor
    matrix6 const P = voigt::kinetic::deviatoric();

private:
    /// Material with micromechanical parameters
    micromechanical_elastomer material;
};

/** \} */
}
