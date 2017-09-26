
#pragma once

#include "HypoElasticPlastic.hpp"

#include "numeric/Tensor.hpp"

#include "material/IsotropicElasticPlastic.hpp"

namespace neon
{
class J2Plasticity : public HypoElasticPlastic
{
public:
    J2Plasticity(InternalVariables& variables, Json::Value const& material_data);

    ~J2Plasticity();

    void update_internal_variables(double const time_step_size) override;

    Material const& intrinsic_material() const override { return material; }

    virtual bool is_finite_deformation() const override { return false; }

protected:
    Matrix3 compute_cauchy_stress(Matrix3 const& elastic_strain) const;

    CMatrix elastic_moduli() const;

    CMatrix deviatoric_projection() const;

    CMatrix algorithmic_tangent(double const plastic_increment,
                                double const accumulated_plastic_strain,
                                double const von_mises,
                                Matrix3 const& n) const;

    /**
     * Performs the radial return algorithm with nonlinear hardening for
     * projecting the stress onto the yield surface.  This provides the plastic
     * increment required for updating the internal variables
     */
    double perform_radial_return(double const von_mises,
                                 double const accumulated_plastic_strain) const;

    /**
     * Evaluates the yield function and returns greater than zero if
     * the yield function has been violated
     */
    double evaluate_yield_function(double const von_mises,
                                   double const accumulated_plastic_strain,
                                   double const plastic_increment = 0.0) const;

protected:
    IsotropicElasticPlastic material;

    CMatrix const C_e = elastic_moduli();
    CMatrix const I_dev = voigt::kinematic::deviatoric();
};
}
