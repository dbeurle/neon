
#pragma once

#include "constitutive/ConstitutiveModel.hpp"

#include <material/LinearElastic.hpp>

// #include <json/forwards.h>

namespace neon::mechanical::solid
{
/**
 * \ingroup Hyperelastic
 * \addtogroup Hyperelastic
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
 *
 */
class compressible_neohooke : public ConstitutiveModel
{
public:
    /**
     * @param variables Reference to internal state variable store
     * @param material_data Json object with material data
     */
    explicit compressible_neohooke(std::shared_ptr<InternalVariables>& variables,
                                   json const& material_data);

    ~compressible_neohooke() = default;

    virtual void update_internal_variables(double const time_step_size) override final;

    virtual Material const& intrinsic_material() const override final { return material; };

    virtual bool is_finite_deformation() const override final { return true; };

private:
    LinearElastic material; //!< Elastic model where C1 = mu/2 and C2 = bulk-modulus / 2
};
/** \} */
}
