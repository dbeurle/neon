
#include "io/json_forward.hpp"

#include <string>

#pragma once

namespace neon
{
class Material
{
public:
    Material(json const& intrinsic_material_data);

    virtual ~Material() = default;

    std::string const& name() const { return material_name; }

    /** @return the density if specified in the material data */
    double initial_density() const;

    /** @return the specific heat if specified in the material data */
    double specific_heat() const;

protected:
    std::string material_name;

    bool is_density_specified = false;
    bool is_specific_heat_specified = false;

    double density_0 = 0.0;
    double c_p = 0.0;
};
}
