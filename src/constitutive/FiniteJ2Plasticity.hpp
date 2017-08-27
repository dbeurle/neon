
#pragma once

#include "J2Plasticity.hpp"

#include "material/IsotropicElasticPlastic.hpp"
#include "numeric/DenseTypes.hpp"
#include "numeric/Tensor.hpp"

#include <array>

namespace neon
{
/**
 * FiniteJ2Plasticity is a material model for the large strain deformation
 * according to the J2 yield theory.  This implementation follows the steps
 * outlined in Computational Methods for Plasticity by Neto et.al, 2008
 * on pages 598.  This does not use the principal stress formulation in order
 * to reuse methods from J2Plasticity.
 *
 * This class is responsible for computation of the Cauchy stress and the material
 * tangent matrix for use in the Newton Raphson iterations.
 *
 * \sa J2Plasticity
 */
class FiniteJ2Plasticity : public J2Plasticity
{
public:
    FiniteJ2Plasticity(InternalVariables& variables, Json::Value const& material_data);

    ~FiniteJ2Plasticity();

    void update_internal_variables(double const time_step_size) override final;

    Material const& intrinsic_material() const override final { return material; }

    virtual bool is_finite_deformation() const override final { return true; };

protected:
    std::tuple<Vector3, std::array<Matrix3, 3>, int> compute_eigenvalues_eigenprojections(
        Matrix3 const& X) const;

    CMatrix derivative_tensor_log(Matrix3 const& Be_trial) const;

    CMatrix derivative_tensor_log_equal(Vector3 const& e,
                                        std::array<Matrix3, 3> const& E,
                                        std::array<int, 3> permutation) const;

    CMatrix deformation_part(Matrix3 const& Be_trial) const;

    CMatrix stress_part(Matrix3 const& cauchy_stress) const;

protected:
    CMatrix const Isym = fourth_order_identity();
};

inline CMatrix FiniteJ2Plasticity::derivative_tensor_log_equal(Vector3 const& x,
                                                               std::array<Matrix3, 3> const& E,
                                                               std::array<int, 3> permutation) const
{
    // BUG Fix the first 2*Isym as placeholder for proper dX^2 / dX expression
    auto const[a, b, c] = permutation;
    return std::log(x(a)) / ((x(a) - x(b)) * (x(a) - x(c)))
               * (2.0 * Isym - (x(b) - x(c)) * Isym - (2.0 * x(a) - x(b) - x(c)) * outer_product(E[a])
                  - (x(b) - x(c)) * (outer_product(E[b]) - outer_product(E[c])))
           + 1.0 / x(a) * outer_product(E[a]);
}
}
