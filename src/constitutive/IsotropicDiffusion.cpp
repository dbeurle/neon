
#include "IsotropicDiffusion.hpp"

#include "InternalVariables.hpp"

namespace neon::diffusion
{
IsotropicDiffusion::IsotropicDiffusion(InternalVariables& variables, Json::Value const& material_data)
    : ConstitutiveModel(variables), material(material_data)
{
    variables.add(InternalVariables::Tensor::Conductivity);

    for (auto& k : variables(InternalVariables::Tensor::Conductivity))
    {
        k = material.conductivity_coefficient() * Matrix3::Identity();
    }
}

void IsotropicDiffusion::update_internal_variables(double const time_step_size) {}
}
