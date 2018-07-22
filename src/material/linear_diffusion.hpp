
#pragma once

#include "material_property.hpp"

namespace neon
{
/// linear_diffusion is responsible for parsing the input file the for the correct
/// material parameters for a thermal simulation.  The heat equation also has
/// additional meaning, for example with specie concentration and modelling
/// diffusive processes.
class linear_diffusion : public material_property
{
public:
    explicit linear_diffusion(json const& material_data);

    [[nodiscard]] auto diffusivity_coefficient() const noexcept
    {
        return conductivity / (specific_heat * density_0);
    }

    [[nodiscard]] auto conductivity_coefficient() const noexcept { return conductivity; }

    [[nodiscard]] auto specific_heat_coefficient() const noexcept { return specific_heat; }

protected:
    double diffusivity{0.0}, conductivity{0.0}, specific_heat{0.0};
};
}
