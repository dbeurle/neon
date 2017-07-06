
#include "Material.hpp"

#include "MaterialExceptions.hpp"

#include <json/value.h>

namespace neon
{
Material::Material(Json::Value const& intrinsic_material_data)
{
    if (intrinsic_material_data["Name"].empty()) throw MaterialPropertyException("Name");

    material_name = intrinsic_material_data["Name"].asString();

    is_density_specified = intrinsic_material_data.isMember("Density");
    if (is_density_specified)
    {
        density_0 = intrinsic_material_data["Density"].asDouble();
    }
}

double Material::initial_density() const
{
    if (!is_density_specified)
        throw std::runtime_error("Density was requested, but not specified in the input file\n");
    return density_0;
}
}
