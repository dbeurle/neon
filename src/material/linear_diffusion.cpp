
#include "linear_diffusion.hpp"

#include "io/json.hpp"

#include "exceptions.hpp"

namespace neon
{
linear_diffusion::linear_diffusion(json const& material_data) : material_property(material_data)
{
    if (material_data.count("Conductivity"))
    {
        conductivity = material_data["Conductivity"];
    }
    if (material_data.count("SpecificHeat"))
    {
        specific_heat = material_data["SpecificHeat"];
    }
    else
    {
        throw std::domain_error("\"Diffusivity\" needs to be specified as a material "
                                "property");
    }
}
}
