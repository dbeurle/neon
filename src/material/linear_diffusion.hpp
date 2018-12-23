
#pragma once

/// @file

#include "material_property.hpp"

namespace neon
{
/// linear_diffusion is responsible for extracting the correct
/// material parameters for a scalar diffusion formulation.
class linear_diffusion : public material_property
{
public:
    explicit linear_diffusion(json const& material_data);

    [[nodiscard]] double diffusivity() const noexcept;

    [[nodiscard]] double conductivity() const noexcept;

    [[nodiscard]] double specific_heat() const noexcept;

protected:
    double m_diffusivity{0.0}, m_conductivity{0.0}, m_specific_heat{0.0};
};
}
