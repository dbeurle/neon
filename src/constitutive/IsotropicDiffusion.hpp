
#pragma once

#include "ConstitutiveModel.hpp"

#include "material/LinearDiffusion.hpp"

namespace neon::diffusion
{
using InternalVariables = neon::InternalVariables<3>;

/**
 * IsotropicDiffusion computes the isotropic constitutive matrix for linear
 * and isotropic diffusion problems.
 */
class IsotropicDiffusion : public diffusion::ConstitutiveModel
{
public:
    IsotropicDiffusion(InternalVariables& variables, Json::Value const& material_data);

    void update_internal_variables(double const time_step_size) override;

    Material const& intrinsic_material() const override { return material; }

    bool is_finite_deformation() const override { return false; }

protected:
    LinearDiffusion material;
};
}
