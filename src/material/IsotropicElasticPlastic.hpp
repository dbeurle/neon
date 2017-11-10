
#pragma once

#include "LinearElastic.hpp"

namespace neon
{
/** PlasticMaterial is a base class for material exhibiting plastic behavior */
class IsotropicElasticPlastic : public LinearElastic
{
public:
    IsotropicElasticPlastic(Json::Value const& material_data);

    ~IsotropicElasticPlastic() = default;

    double yield_stress(double const effective_strain) const;

    double hardening_modulus(double const effective_strain) const;

    double kinematic_modulus(double const effective_strain) const;

protected:
    double stress_y = 0.0;

    double H_iso = 0.0; //!< Isotropic hardening modulus
    double K_iso = 0.0; //!< Kinematic hardening modulus
};
}
