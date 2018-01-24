
#pragma once

#include "IsotropicElasticPlastic.hpp"

namespace neon
{
/** PlasticMaterial is a base class for material exhibiting plastic behavior */
class IsotropicElasticPlasticDamage : public IsotropicElasticPlastic
{
public:
    IsotropicElasticPlasticDamage(json const& material_data);

    ~IsotropicElasticPlasticDamage() = default;

    double softening_multiplier() const { return gamma; }
    double kinematic_hardening_modulus() const { return C; }
    double plasticity_viscous_exponent() const { return np; }
    double plasticity_viscous_multiplier() const { return kp; }
    double damage_viscous_exponent() const { return nd; }
    double damage_viscous_multiplier() const { return kd; }

protected:
    double gamma = 1.0; // !< Kinematic hardening numerator
    double C = 1.0;     // !< Kinematic hardening denominator
    double kp = 1.0;    //!< viscous multiplier
    double np = 1.0;    //!< viscous exponent
    double kd = 1.0;    //!< damage viscous multiplier
    double nd = 1.0;    //!< damage viscous exponent
};
}
