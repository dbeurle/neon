
#pragma once

#include "small_strain_J2_plasticity.hpp"

#include "numeric/dense_matrix.hpp"
#include "numeric/tensor_operations.hpp"

#include <array>

namespace neon::mechanical::solid
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
    finite_strain_J2_plasticity(std::shared_ptr<InternalVariables>& variables,
                                json const& material_data);

    ~finite_strain_J2_plasticity();

    void update_internal_variables(double const time_step_size) override final;

    Material const& intrinsic_material() const override final { return material; }

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
     * @param J Determinant of the deformation gradient
     * @param Be_trial Elastic trial left Cauchy Green deformation tensor
     * @param cauchy_stress Cauchy stress
     * @param D Tangent matrix from small-strain theory
     *
     * \sa compute_B
     * \sa compute_L
     * \sa compute_H
     */
    matrix6 consistent_tangent(double const J,
                               matrix3 const& Be_trial,
                               matrix3 const& cauchy_stress,
                               matrix6 const& D);

    /**
     * Computes the fourth order B matrix
     * \f{align*}{
        B_{ijkl} &= \delta_{ik} (\mathbf{B}_{n+1}^{e, trial})_{jl} + \delta_{jk}
     (\mathbf{B}_{n+1}^{e, trial})_{il} \f}
     */
    matrix6 compute_B(matrix3 const& Be_trial) const;

    /**
     * Computes the derivative of the tensor log with respect to the elastic trial
     * left Cauchy-Green tensor
       \f{align*}{
         L &= \frac{\partial \ln \mathbf{B}_e^{trial}}{\partial \mathbf{B}_e^{trial}}
       \f}
     */
    matrix6 compute_L(matrix3 const& Be_trial) const;

    /**
     * Computes the finite strain geometric contribution
       \f{align*}{
        H_{ijkl} &= \sigma_{il}\delta_{jk} + \sigma_{jl}\delta_{ik}
       \f}
     */
    matrix6 compute_H(matrix3 const& cauchy_stress) const;

    /**
     *
     */
    std::tuple<vector3, std::array<matrix3, 3>, bool, std::array<int, 3>> compute_eigenvalues_eigenprojections(
        matrix3 const& X) const;

    /**
     * Computes the derivative when all the eigenvalues are unique
     */
    matrix6 derivative_tensor_log_unique(matrix3 const& Be_trial,
                                         vector3 const& e,
                                         std::array<matrix3, 3> const& E,
                                         std::array<int, 3> const& abc_ordering) const;

    /**
     * Computes the derivative of the tensor squared with respect to the
     * same tensor
       \f{align*}{
        \mathbf{D}(\mathbf{X}) &= \frac{\partial \mathbf{X}^2}{\partial \mathbf{X}}
       \f}
     */
    matrix6 dX2_dX(matrix3 const& X) const;

    /**
     * Transforms form a fourth-order tensor in Voigt notation
     * to a fourth-order tensor in mandel notation
     * \sa mandel_transformation
     */
    matrix6 voigt_to_mandel(matrix6 const& V) const { return V.array() * M.array(); }

    /**
     * Provides the transformation between Voigt notation and Mandel notation
     * for providing a method to perform double dot tensor products between
     * fourth order tensors
     */
    matrix6 mandel_transformation() const;

protected:
    matrix6 const Isym = voigt::kinematic::fourth_order_identity();
    matrix6 const M = mandel_transformation();
};
}
