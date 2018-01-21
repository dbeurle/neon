
#pragma once

#include "constitutive/ConstitutiveModel.hpp"

#include <material/LinearElastic.hpp>

#include <json/forwards.h>

namespace neon::mechanical::solid
{
class NeoHooke : public ConstitutiveModel
{
public:
    /**
     * @param variables Reference to internal state variable store
     * @param material_data Json object with material data
     */
    explicit NeoHooke(std::shared_ptr<InternalVariables>& variables,
                      json const& material_data);

    ~NeoHooke() = default;

    virtual void update_internal_variables(double const time_step_size) override final;

    virtual Material const& intrinsic_material() const override final { return material; };

    virtual bool is_finite_deformation() const override final { return true; };

private:
    LinearElastic material; //!< Elastic model where C1 = mu/2 and C2 = bulk-modulus / 2
};
}
