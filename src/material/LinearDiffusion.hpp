
#pragma once

#include "Material.hpp"

#include <json/forwards.h>

namespace neon
{
class LinearDiffusion : public Material
{
public:
    LinearDiffusion(Json::Value const& material_data);

    auto diffusivity_coefficient() const { return conductivity / (Cp * density_0); }

    auto conductivity_coefficient() const { return conductivity; }

    auto specific_heat() const { return Cp; }

protected:
    double diffusivity = 0.0;
    double conductivity = 0.0;
    double Cp = 0.0; // Specific heat
};
}
