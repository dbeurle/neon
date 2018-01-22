
#include "LinearDiffusion.hpp"

#include "io/json.hpp"

#include "Exceptions.hpp"

namespace neon
{
LinearDiffusion::LinearDiffusion(json const& material_data) : Material(material_data)
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
        throw MaterialPropertyException("\"Diffusivity\" needs to be specified as a material "
                                        "property");
    }
}
}
