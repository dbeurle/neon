
#include "LinearDiffusion.hpp"

#include <json/json.h>

#include "Exceptions.hpp"

namespace neon
{
LinearDiffusion::LinearDiffusion(Json::Value const& material_data) : Material(material_data)
{
    if (material_data.isMember("Conductivity"))
    {
        conductivity = material_data["Conductivity"].asDouble();
    }
    if (material_data.isMember("SpecificHeat"))
    {
        Cp = material_data["SpecificHeat"].asDouble();
    }
    else
    {
        throw MaterialPropertyException("\"Diffusivity\" needs to be specified as a material "
                                        "property");
    }
}
}
