
#pragma once

#include "J2Plasticity.hpp"

#include "material/IsotropicElasticPlastic.hpp"
#include "numeric/DenseTypes.hpp"
#include "numeric/Tensor.hpp"

#include <array>

namespace neon
{
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

    CMatrix deformation_part(Matrix3 const& Be_trial) const;

    CMatrix stress_part(Matrix3 const& cauchy_stress) const;

protected:
    CMatrix const Isym = fourth_order_identity();
};
}
