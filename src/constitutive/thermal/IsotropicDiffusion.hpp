
#pragma once

#include "constitutive/ConstitutiveModel.hpp"

#include "material/LinearDiffusion.hpp"

namespace neon::diffusion
{
/**
 * IsotropicDiffusion computes the isotropic constitutive matrix for linear
 * and isotropic diffusion problems.
 */
class IsotropicDiffusion : public diffusion::ConstitutiveModel
{
public:
    IsotropicDiffusion(std::shared_ptr<InternalVariables>& variables,
                       Json::Value const& material_data);

    void update_internal_variables(double const time_step_size) override;

    Material const& intrinsic_material() const override { return material; }

    bool is_finite_deformation() const override { return false; }

protected:
    LinearDiffusion material;
};
}
