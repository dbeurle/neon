
#pragma once

#include "constitutive/ConstitutiveModel.hpp"

#include "numeric/DenseMatrix.hpp"

#include "material/LinearElastic.hpp"

namespace neon::mechanical::solid
{
/**
 * IsotropicLinearElasticity is responsible for compute the moduli and the
 * stress for the three-dimensional theory.  See \cite Hughes2012 for the
 * theoretical developments.
 */
class IsotropicLinearElasticity : public ConstitutiveModel
{
public:
    explicit IsotropicLinearElasticity(std::shared_ptr<InternalVariables>& variables,
                                       Json::Value const& material_data);

    virtual ~IsotropicLinearElasticity();

    virtual void update_internal_variables(double const time_step_size) override;

    [[nodiscard]] virtual Material const& intrinsic_material() const override { return material; }

    [[nodiscard]] virtual bool is_finite_deformation() const override { return false; }

protected:
    [[nodiscard]] Matrix6 elastic_moduli() const;

private:
    LinearElastic material;

protected:
    Matrix6 const C_e = elastic_moduli();
};
}
