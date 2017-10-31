
#pragma once

#include "Material.hpp"

#include <json/forwards.h>

namespace neon
{
class LinearElastic : public Material
{
public:
    explicit LinearElastic(Json::Value const& material_data);

    /** @return Elastic modulus or Young's modulus */
    auto elastic_modulus() const { return 9.0 * K * G / (3.0 * K + G); }

    /** @return The Poisson's ratio */
    auto Poissons_ratio() const { return (3.0 * K - 2.0 * G) / (2.0 * (3.0 * K + G)); }

    /** Compute the Lame parameter \f[ \mu = \frac{E}{2(1+\nu)} \f] */
    auto mu() const { return G; }

    /** Compute the Lame parameter \f[ \lambda = \frac{\nu E}{(1+\nu)(1-2\nu)} \f] */
    auto lambda() const { return K - 2.0 * G / 3.0; }

    /** @return Bulk modulus \sa lambda \sa mu */
    auto bulk_modulus() const { return K; }

    /** @return Shear modulus \sa mu */
    auto shear_modulus() const { return G; }

    /** @return a pair of Lame's parameters lambda and mu respectively */
    auto Lame_parameters() const { return std::make_pair(lambda(), mu()); }

protected:
    double K; //!< Bulk modulus
    double G; //!< Shear modulus
};
}
