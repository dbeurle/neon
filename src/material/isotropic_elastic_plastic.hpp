
#pragma once

#include "isotropic_elastic_property.hpp"

namespace neon
{
/** PlasticMaterial is a base class for material exhibiting plastic behavior */
class isotropic_elastic_plastic : public isotropic_elastic_property
{
public:
    explicit isotropic_elastic_plastic(json const& material_data);

    ~isotropic_elastic_plastic() = default;

    [[nodiscard]] double yield_stress(double const effective_strain) const;

    [[nodiscard]] double hardening_modulus(double) const;

    [[nodiscard]] double kinematic_modulus(double) const;

protected:
    /// Yield stress
    double stress_y = 0.0;
    /// Isotropic hardening modulus
    double H_iso = 0.0;
    /// Kinematic hardening modulus
    double K_iso = 0.0;
};
}
