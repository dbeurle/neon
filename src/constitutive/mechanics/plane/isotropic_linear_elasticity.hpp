
#pragma once

#include "constitutive/constitutive_model.hpp"

#include "numeric/dense_matrix.hpp"

#include "material/isotropic_elastic_property.hpp"

namespace neon::mechanics::plane
{
class isotropic_linear_elasticity : public constitutive_model
{
public:
    enum class plane { stress, strain };

public:
    /** Provide an internal variable class to be populated by the constitutive model */
    explicit isotropic_linear_elasticity(std::shared_ptr<internal_variables_t>& variables,
                                         json const& material_data,
                                         plane const state);

    ~isotropic_linear_elasticity();

    /**
     * Update the required internal variables and tangent matrix at quadrature
     * points
     * @param time_step_size Time step size (or load increment if quasi-static)
     */
    virtual void update_internal_variables(double const time_step_size);

    /** @return A base class reference to the common material properties */
    [[nodiscard]] virtual material_property const& intrinsic_material() const { return material; }

    [[nodiscard]] virtual bool is_finite_deformation() const { return false; };

    [[nodiscard]] virtual bool is_symmetric() const { return true; };

protected:
    [[nodiscard]] matrix3 elastic_moduli() const;

private:
    isotropic_elastic_property material;

protected:
    plane state;

    matrix3 const C_e = elastic_moduli();
};
}
