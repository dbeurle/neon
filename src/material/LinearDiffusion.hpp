
#pragma once

#include "Material.hpp"

#include <json/forwards.h>

namespace neon
{
/**
 * LinearDiffusion is responsible for parsing the input file the for the correct
 * material parameters for a thermal simulation.  The heat equation also has
 * additional meaning, for example with specie concentration and modelling
 * diffusive processes.
 */
class LinearDiffusion : public Material
{
public:
    LinearDiffusion(json const& material_data);

    auto diffusivity_coefficient() const { return conductivity / (specific_heat * density_0); }

    auto conductivity_coefficient() const { return conductivity; }

    auto specific_heat_coefficient() const { return specific_heat; }

protected:
    double diffusivity = 0.0, conductivity = 0.0, specific_heat = 0.0;
};
}
