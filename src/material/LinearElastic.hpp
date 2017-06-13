
#pragma once

#include "Material.hpp"

#include <json/forwards.h>

namespace neon
{
class LinearElastic : public Material
{
public:
    LinearElastic(Json::Value const& material_data);

    /** @return The elastic modulus */
    auto elastic_modulus() const { return E; }

    /** @return The Poisson's ratio */
    auto Poissons_ratio() const { return nu; }

    /** Compute the Lame parameter \f[ \mu = \frac{E}{2(1+\nu)} \f] */
    auto mu() const { return E / (2.0 * (1.0 + nu)); }

    /** Compute the Lame parameter \f[ \lambda = \frac{\nu E}{(1+\nu)(1-2\nu)} \f] */
    auto lambda() const { return E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu)); }

    auto bulk_modulus() const { return lambda() + 2.0 / 3.0 * mu(); }

    /** @return A pair with mu and lambda respectively */
    auto LameConstants() const { return std::make_pair(mu(), lambda()); }

protected:
    double E;  //!< Elastic modulus
    double nu; //!< Poisson's ratio
};
}
