
#pragma once

/// @file

#include "isotropic_linear_elasticity.hpp"

#include "numeric/tensor_operations.hpp"

#include "material/isotropic_elastic_plastic.hpp"

namespace neon::mechanics::solid
{
/**
 * small_strain_J2_plasticity is responsible for computing the small strain J2 plasticity
 * stress and tangent operator matrix for each internal variable at the quadrature
 * points.  The method can be found for small strain plasticity using a combination
 * of @cite Neto2011 and @cite Belytschko2013nonlinear
 */
class small_strain_J2_plasticity : public isotropic_linear_elasticity
{
public:
    explicit small_strain_J2_plasticity(std::shared_ptr<internal_variables_t>& variables,
                                        json const& material_data);

    virtual ~small_strain_J2_plasticity();

    virtual void update_internal_variables(double) override;

    virtual material_property const& intrinsic_material() const override { return material; }

    virtual bool is_finite_deformation() const override { return false; }

protected:
    [[nodiscard]] matrix6 algorithmic_tangent(double const plastic_increment,
                                              double const accumulated_plastic_strain,
                                              double const von_mises,
                                              matrix3 const& normal) const;

    /**
     * Performs the radial return algorithm with nonlinear hardening for
     * projecting the stress onto the yield surface.  This provides the plastic
     * increment required for updating the internal variables
     */
    [[nodiscard]] double perform_radial_return(double const von_mises,
                                               double const accumulated_plastic_strain) const;

    /**
     * Evaluates the yield function and returns greater than zero if
     * the yield function has been violated
     */
    [[nodiscard]] double evaluate_yield_function(double const von_mises,
                                                 double const accumulated_plastic_strain,
                                                 double const plastic_increment = 0.0) const;

protected:
    isotropic_elastic_plastic material;

    matrix6 const I_dev = voigt::kinematic::deviatoric();
};
}
