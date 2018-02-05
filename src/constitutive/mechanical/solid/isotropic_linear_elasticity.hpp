
#pragma once

#include "constitutive/ConstitutiveModel.hpp"

#include "numeric/dense_matrix.hpp"

#include "material/LinearElastic.hpp"

namespace neon::mechanical::solid
{
/**
 * isotropic_linear_elasticity is responsible for compute the moduli and the
 * stress for the three-dimensional theory.  See \cite Hughes2012 for the
 * theoretical developments.
 */
class isotropic_linear_elasticity : public ConstitutiveModel
{
public:
    explicit isotropic_linear_elasticity(std::shared_ptr<InternalVariables>& variables,
                                       json const& material_data);

    virtual ~isotropic_linear_elasticity();

    virtual void update_internal_variables(double const time_step_size) override;

    [[nodiscard]] virtual Material const& intrinsic_material() const override { return material; }

    [[nodiscard]] virtual bool is_finite_deformation() const override { return false; }

protected:
    [[nodiscard]] matrix6 elastic_moduli() const;

private:
    LinearElastic material;

protected:
    matrix6 const C_e = elastic_moduli();
};
}
