
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

    double yield_stress(double effective_strain) const;

    double hardening_modulus(double effective_strain) const;

    double kinematic_modulus(double effective_strain) const;

protected:
    double stress_y = 0.0;

    double H = 0.0; //!< Isotropic hardening modulus
    double K = 0.0; //!< Kinematic hardening modulus
};
}
