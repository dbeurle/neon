
#pragma once

namespace neon
{
class InternalVariables;
class PlasticMaterial;

/**
 * ConstitutiveModel is the base class for all constitutive models.  The derived
 * classes define their own internal variables, implement an update internal
 * variables routine and update a constitutive model for use in the global
 * assembly routine
 */
class ConstitutiveModel
{
public:
    /** Provide an internal variable class to be populated by the constitutive model */
    ConstitutiveModel(InternalVariables& variables) : variables(variables) {}

    /** Update the required internal variables at quadrature points */
    virtual void update_internal_variables() = 0;

    /** Update constitutive matrix at quadrature points */
    virtual void update_continuum_tangent() = 0;

protected:
    InternalVariables& variables;
};
}
