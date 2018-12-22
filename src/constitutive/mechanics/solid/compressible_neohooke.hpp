
#pragma once

/// @file

#include "constitutive/constitutive_model.hpp"

#include <material/isotropic_elastic_property.hpp>

namespace neon::mechanics::solid
{
/**
 * \ingroup Hyperelastic
 * \{
 *
 * The compressible_neohooke material is a hyperelastic material model that
 * uses a volumetric free energy function in addition to the micromechanical
 * contribution
   \f{align*}{
       \psi &= \frac{1}{2}\mu_0 (I_\boldsymbol{C} - 3)
             + \frac{1}{2}\lambda_0 (\ln{J})^2 - \mu_0 \ln J
   \f}
 * where \f$ \lambda_0 \f$ and \f$ \mu \f$ are the Lamé parameters.
 */
class compressible_neohooke : public constitutive_model
{
public:
    /// \param variables Reference to internal state variable store
    /// \param material_data Json object with material data
    explicit compressible_neohooke(std::shared_ptr<internal_variables_t>& variables,
                                   json const& material_data);

    virtual ~compressible_neohooke() = default;

    void update_internal_variables(double) override final;

    material_property const& intrinsic_material() const override final { return material; };

    bool is_finite_deformation() const override final { return true; };

private:
    /// Elastic model where C1 = mu/2 and C2 = bulk-modulus / 2
    isotropic_elastic_property material;
};
/** \} */
}
