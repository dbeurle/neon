
#pragma once

/// @file

#include "constitutive/constitutive_model.hpp"
#include "material/linear_diffusion.hpp"

namespace neon::diffusion
{
/// isotropic_diffusion computes the isotropic constitutive matrix for linear
/// and isotropic diffusion problems.
class isotropic_diffusion : public constitutive_model
{
public:
    isotropic_diffusion(std::shared_ptr<internal_variables_t>& variables, json const& material_data);

    void update_internal_variables(double const time_step_size) override;

    material_property const& intrinsic_material() const override { return material; }

    bool is_finite_deformation() const override final { return false; }

protected:
    linear_diffusion material;
};
}
