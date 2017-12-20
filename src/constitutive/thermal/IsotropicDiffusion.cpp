
#include "IsotropicDiffusion.hpp"

#include "constitutive/InternalVariables.hpp"

namespace neon::diffusion
{
IsotropicDiffusion::IsotropicDiffusion(std::shared_ptr<InternalVariables>& variables,
                                       Json::Value const& material_data)
    : ConstitutiveModel(variables), material(material_data)
{
    variables->add(InternalVariables::Tensor::Conductivity);

    for (auto& k : variables->fetch(InternalVariables::Tensor::Conductivity))
    {
        k = material.conductivity_coefficient() * matrix3::Identity();
    }
}

void IsotropicDiffusion::update_internal_variables(double const time_step_size) {}
}
