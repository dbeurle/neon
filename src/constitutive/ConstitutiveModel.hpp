
#pragma once

#include "InternalVariablesForwards.hpp"

#include <memory>

namespace neon
{
class material_property;

/**
 * ConstitutiveModel is the templated base class for all constitutive models.
 * The derived classes define their own internal variables, implement an update
 * internal variables routine and update a constitutive model for use in the global
 * assembly routine
 */
template <int rank2_dimension, int rank4_dimension>
class ConstitutiveModel
{
public:
    /** Provide an internal variable class to be populated by the constitutive model */
    explicit ConstitutiveModel(
        std::shared_ptr<InternalVariables<rank2_dimension, rank4_dimension>>& variables)
        : variables(variables)
    {
    }

    /**
     * Update the required internal variables and tangent matrix at quadrature
     * points
     * @param time_step_size Time step size (or load increment if quasi-static)
     */
    virtual void update_internal_variables(double const time_step_size) = 0;

    /** @return A base class reference to the common material properties */
    [[nodiscard]] virtual material_property const& intrinsic_material() const = 0;

    [[nodiscard]] virtual bool is_finite_deformation() const = 0;

    [[nodiscard]] virtual bool is_symmetric() const { return true; };

protected:
    std::shared_ptr<InternalVariables<rank2_dimension, rank4_dimension>> variables;
};

namespace mechanical::solid
{
using ConstitutiveModel = neon::ConstitutiveModel<3, 6>;
}
namespace mechanical::plane
{
using ConstitutiveModel = neon::ConstitutiveModel<2, 3>;
}
namespace diffusion
{
using ConstitutiveModel = neon::ConstitutiveModel<3, 3>;
}
}
