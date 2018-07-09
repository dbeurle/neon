
#pragma once

#include "material/material_property.hpp"

namespace neon
{
class isotropic_elastic_property : public material_property
{
public:
    explicit isotropic_elastic_property(json const& material_data);

    /// \return Elastic modulus or Young's modulus
    auto elastic_modulus() const noexcept { return 9.0 * K * G / (3.0 * K + G); }

    /// \return Poisson's ratio
    auto Poissons_ratio() const noexcept { return (3.0 * K - 2.0 * G) / (2.0 * (3.0 * K + G)); }

    /// Compute the Lame parameter \f[ \mu = \frac{E}{2(1+\nu)} \f]
    auto mu() const noexcept { return G; }

    /// Compute the Lame parameter \f[ \lambda = \frac{\nu E}{(1+\nu)(1-2\nu)} \f]
    auto lambda() const noexcept { return K - 2.0 * G / 3.0; }

    /// \return Bulk modulus \sa lambda \sa mu
    auto bulk_modulus() const noexcept { return K; }

    /// \return Shear modulus \sa mu
    auto shear_modulus() const noexcept { return G; }

    /// \return a pair of Lame's parameters lambda and mu respectively
    auto Lame_parameters() const noexcept { return std::make_pair(lambda(), mu()); }

protected:
    /// Bulk modulus
    double K;
    /// Shear modulus
    double G;
};
}
