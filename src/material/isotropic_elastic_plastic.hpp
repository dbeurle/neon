
#pragma once

#include "isotropic_elastic_property.hpp"

namespace neon
{
/** PlasticMaterial is a base class for material exhibiting plastic behavior */
class isotropic_elastic_plastic : public isotropic_elastic_property
{
public:
    isotropic_elastic_plastic(json const& material_data);

    ~isotropic_elastic_plastic() = default;

    double yield_stress(double const effective_strain) const;

    double hardening_modulus(double const effective_strain) const;

    double kinematic_modulus(double const effective_strain) const;

protected:
    double stress_y = 0.0;

    double H_iso = 0.0; //!< Isotropic hardening modulus
    double K_iso = 0.0; //!< Kinematic hardening modulus
};
}
