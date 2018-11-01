
#include "material_property.hpp"

#include <stdexcept>

#include "io/json.hpp"

namespace neon
{
material_property::material_property(json const& intrinsic_material_data)
{
    if (!intrinsic_material_data.count("name"))
    {
        throw std::domain_error("\"name\" must be specified for a material");
    }

    material_name = intrinsic_material_data["name"];

    is_density_specified = intrinsic_material_data.count("density");

    if (is_density_specified)
    {
        density_0 = intrinsic_material_data["density"];
    }

    is_specific_heat_specified = intrinsic_material_data.count("specific_heat");

    if (is_specific_heat_specified)
    {
        c_p = intrinsic_material_data["specific_heat"];
    }
}

double material_property::initial_density() const
{
    if (!is_density_specified)
    {
        throw std::domain_error("density was requested, but not specified in the input "
                                "file\n");
    }
    return density_0;
}

double material_property::specific_heat() const
{
    if (!is_specific_heat_specified)
    {
        throw std::domain_error("specific_heat was requested, but not specified in the input "
                                "file\n");
    }
    return c_p;
}
}
