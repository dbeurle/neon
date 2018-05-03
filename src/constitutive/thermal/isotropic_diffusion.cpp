
#include "isotropic_diffusion.hpp"
#include "constitutive/internal_variables.hpp"

namespace neon::diffusion
{
isotropic_diffusion::isotropic_diffusion(std::shared_ptr<internal_variables_t>& variables,
                                         json const& material_data)
    : constitutive_model(variables), material(material_data)
{
    variables->add(internal_variables_t::second::Conductivity);

    for (auto& k : variables->fetch(internal_variables_t::second::Conductivity))
    {
        k = material.conductivity_coefficient() * matrix3::Identity();
    }
}

void isotropic_diffusion::update_internal_variables(double const time_step_size) {}
}
