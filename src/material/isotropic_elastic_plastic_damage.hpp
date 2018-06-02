
#pragma once

#include "isotropic_elastic_plastic.hpp"

namespace neon
{
/**
 * isotropic_elastic_plastic_damage is responsible for providing an interface
 * for the parameters required in an isotropic damage model.
 * \sa small_strain_J2_plasticity_damage
 */
class isotropic_elastic_plastic_damage : public isotropic_elastic_plastic
{
public:
    isotropic_elastic_plastic_damage(json const& material_data);

    double softening_multiplier() const { return gamma; }
    double kinematic_hardening_modulus() const { return C; }
    double plasticity_viscous_exponent() const { return np; }
    double plasticity_viscous_denominator() const { return sp; }
    double damage_viscous_exponent() const { return nd; }
    double damage_viscous_denominator() const { return sd; }

protected:
    double gamma = 1.0; // !< Kinematic hardening numerator
    double C = 1.0;     // !< Kinematic hardening denominator
    double sp = 1.0;    //!< viscous denominator
    double np = 1.0;    //!< viscous exponent
    double sd = 1.0;    //!< damage viscous denominator
    double nd = 1.0;    //!< damage viscous exponent
};
}
