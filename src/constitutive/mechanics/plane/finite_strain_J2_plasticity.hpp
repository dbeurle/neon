
#pragma once

#include "constitutive/mechanics/plane/small_strain_J2_plasticity.hpp"

#include "material/isotropic_elastic_plastic.hpp"
#include "numeric/dense_matrix.hpp"
#include "numeric/tensor_operations.hpp"

namespace neon::mechanics::plane
{
/**
 * finite_strain_J2_plasticity is a material model for the large strain deformation
 * according to the J2 yield theory.  This implementation follows the steps
 * outlined in \cite Neto2011 on pages 598.  This does not use the principal
 * stress formulation in order to reuse methods from small_strain_J2_plasticity.
 *
 * This class is responsible for computation of the Cauchy stress and the material
 * tangent matrix for use in the Newton Raphson iterations.
 *
 * \sa small_strain_J2_plasticity
 */
class finite_strain_J2_plasticity : public small_strain_J2_plasticity
{
public:
    finite_strain_J2_plasticity(std::shared_ptr<internal_variables_t>& variables,
                                json const& material_data);

    ~finite_strain_J2_plasticity();

    void update_internal_variables(double) override final;

    material_property const& intrinsic_material() const override final { return material; }

    virtual bool is_finite_deformation() const override final { return true; };

protected:
    /**
     * Computes the tangent modulus \f$\mathbf{c}\f$ for use in the small strain material tangent
     * computation routines
       \f{align*}{
        c_{ijkl} &= a_{ijkl} − \sigma_{jl} \delta_{ik} \;\;\; \text{where} \\
        a_{ijkl} &= \frac{1}{2J} (\mathbf{D} : \mathbf{L} : \mathbf{B})_{ijkl} -
       \sigma_{il}\delta_{jk}
       \f}
       This can be simplified to
       \f{align*}{
          c_{ijkl} &=\frac{1}{2J} (\mathbf{D} : \mathbf{L} : \mathbf{B})_{ijkl} - H_{ijkl}
       \f}
     * where \f$ H_{ijkl} = \sigma_{il}\delta_{jk} + \sigma_{jl}\delta_{ik} \f$ for the
     * Newton-Raphson iterations.
     *
     *
     * \param J Determinant of the deformation gradient
     * \param Be_trial Elastic trial left Cauchy Green deformation tensor
     * \param cauchy_stress Cauchy stress
     * \param D Tangent matrix from small-strain theory
     *
     * \return Consistent tangent operator
     *
     * \sa compute_B
     * \sa compute_L
     * \sa compute_H
     */
    matrix3 consistent_tangent(double const J,
                               matrix2 const& Be_trial,
                               matrix2 const& cauchy_stress,
                               matrix3 const& tangent_operator);

    /**
     * Computes the fourth order B matrix
     * \f{align*}{
        B_{ijkl} &= \delta_{ik} (\mathbf{B}_{n+1}^{e, trial})_{jl} + \delta_{jk}
     (\mathbf{B}_{n+1}^{e, trial})_{il} \f}
     */
    matrix3 compute_B(matrix2 const& Be_trial) const;

    /**
     * Computes the derivative of the tensor log with respect to the elastic trial
     * left Cauchy-Green tensor
       \f{align*}{
         L &= \frac{\partial \ln \mathbf{B}_e^{trial}}{\partial \mathbf{B}_e^{trial}}
       \f}
     * where \f$ \mathbf{B}_e^{trial} \f$ is the elastic trial left Cauchy
     * Green deformation tensor
     */
    matrix3 compute_L(matrix2 const& trial_elastic_left_cauchy_green) const;

    /**
     * Computes the finite strain geometric contribution
       \f{align*}{
        H_{ijkl} &= \sigma_{il}\delta_{jk} + \sigma_{jl}\delta_{ik}
       \f}
     */
    matrix3 compute_H(matrix2 const& cauchy_stress) const;

protected:
    matrix3 const Isym = voigt::kinematic::d2::fourth_order_identity();
};
}
