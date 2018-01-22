
#include "Material.hpp"

#include "Exceptions.hpp"

#include "io/json.hpp"

namespace neon
{
Material::Material(json const& intrinsic_material_data)
{
    if (!intrinsic_material_data.count("Name"))
    {
        throw MaterialPropertyException("Name");
    }

    material_name = intrinsic_material_data["Name"];

    is_density_specified = intrinsic_material_data.count("Density");

    if (is_density_specified)
    {
        density_0 = intrinsic_material_data["Density"];
    }

    is_specific_heat_specified = intrinsic_material_data.count("SpecificHeat");

    if (is_specific_heat_specified)
    {
        c_p = intrinsic_material_data["SpecificHeat"];
    }
}

double Material::initial_density() const
{
    if (!is_density_specified)
    {
        throw std::runtime_error("Density was requested, but not specified in the input "
                                 "file\n");
    }
    return density_0;
}

double Material::specific_heat() const
{
    if (!is_specific_heat_specified)
    {
        throw std::runtime_error("Density was requested, but not specified in the input "
                                 "file\n");
    }
    return c_p;
}
}
