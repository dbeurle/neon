
#pragma once

#include "Material.hpp"

#include <json/forwards.h>

namespace neon
{
class LinearDiffusion : public Material
{
public:
    LinearDiffusion(Json::Value const& material_data);

    auto diffusivity_coefficient() const { return diffusivity; }

protected:
    double diffusivity = 0.0;
};
}
