
#pragma once

#include "constitutive/ConstitutiveModel.hpp"

#include "numeric/DenseMatrix.hpp"

#include "material/LinearElastic.hpp"

namespace neon::mechanical::plane
{
class IsotropicLinearElasticity : public ConstitutiveModel
{
public:
    enum class State { PlaneStress, PlaneStrain };

public:
    /** Provide an internal variable class to be populated by the constitutive model */
    IsotropicLinearElasticity(InternalVariables& variables,
                              Json::Value const& material_data,
                              State const state);

    ~IsotropicLinearElasticity();

    /**
     * Update the required internal variables and tangent matrix at quadrature
     * points
     * @param time_step_size Time step size (or load increment if quasi-static)
     */
    virtual void update_internal_variables(double const time_step_size);

    /** @return A base class reference to the common material properties */
    [[nodiscard]] virtual Material const& intrinsic_material() const { return material; }

    [[nodiscard]] virtual bool is_finite_deformation() const { return false; };

    [[nodiscard]] virtual bool is_symmetric() const { return true; };

protected:
    [[nodiscard]] Matrix2 compute_cauchy_stress(Matrix2 const& elastic_strain) const;

    [[nodiscard]] Matrix3 elastic_moduli() const;

private:
    LinearElastic material;

protected:
    State state;

    Matrix3 const C_e = elastic_moduli();
};
}
