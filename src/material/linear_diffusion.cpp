
#include "linear_diffusion.hpp"

#include "io/json.hpp"

#include "exceptions.hpp"

namespace neon
{
linear_diffusion::linear_diffusion(json const& material_data) : material_property(material_data)
{
    if (material_data.find("conductivity") != material_data.end())
    {
        m_conductivity = material_data["conductivity"];
    }
    if (material_data.find("specific_heat") != material_data.end())
    {
        m_specific_heat = material_data["specific_heat"];
    }
    else
    {
        throw std::domain_error("\"specific_heat\" needs to be specified as a material "
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
