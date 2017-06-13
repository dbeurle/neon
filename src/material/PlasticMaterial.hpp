
#pragma once

#include "LinearElastic.hpp"

namespace neon
{
/** PlasticMaterial is a base class for material exhibiting plastic behavior */
class PlasticMaterial : public LinearElastic
{
public:
    PlasticMaterial(Json::Value const& material_data) : LinearElastic(material_data) {}

    virtual ~PlasticMaterial() = default;

    virtual double evaluate_yield_function(double effective_strain) const = 0;

    virtual double plastic_modulus(double effective_strain) const = 0;
};
}
