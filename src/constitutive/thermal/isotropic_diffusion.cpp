
#include "isotropic_diffusion.hpp"

#include "constitutive/InternalVariables.hpp"

namespace neon::diffusion
{
isotropic_diffusion::isotropic_diffusion(std::shared_ptr<InternalVariables>& variables,
                                         json const& material_data)
    : ConstitutiveModel(variables), material(material_data)
{
    variables->add(InternalVariables::Tensor::Conductivity);

    for (auto& k : variables->fetch(InternalVariables::Tensor::Conductivity))
    {
        k = material.conductivity_coefficient() * matrix3::Identity();
    }
}

void isotropic_diffusion::update_internal_variables(double const time_step_size) {}
}
