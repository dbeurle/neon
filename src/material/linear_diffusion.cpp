
#include "linear_diffusion.hpp"

#include "io/json.hpp"

#include "exceptions.hpp"

namespace neon
{
linear_diffusion::linear_diffusion(json const& material_data) : material_property(material_data)
{
    if (material_data.find("Conductivity") != material_data.end())
    {
        m_conductivity = material_data["Conductivity"];
    }
    if (material_data.find("SpecificHeat") != material_data.end())
    {
        m_specific_heat = material_data["SpecificHeat"];
    }
    else
    {
        throw std::domain_error("\"Diffusivity\" needs to be specified as a material "
                                "property");
    }
}

double linear_diffusion::diffusivity() const noexcept
{
    return m_conductivity / (m_specific_heat * density_0);
}

double linear_diffusion::conductivity() const noexcept { return m_conductivity; }

double linear_diffusion::specific_heat() const noexcept { return m_specific_heat; }
}
