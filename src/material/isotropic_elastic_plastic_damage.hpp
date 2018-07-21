
#pragma once

#include "isotropic_elastic_plastic.hpp"

namespace neon
{
/// isotropic_elastic_plastic_damage is responsible for providing an interface
/// for the parameters required in an isotropic damage model.
/// \sa small_strain_J2_plasticity_damage
class isotropic_elastic_plastic_damage : public isotropic_elastic_plastic
{
public:
    explicit isotropic_elastic_plastic_damage(json const& material_data);

    [[nodiscard]] double softening_multiplier() const noexcept { return gamma; }

    [[nodiscard]] double kinematic_hardening_modulus() const noexcept { return C; }

    [[nodiscard]] double plasticity_viscous_exponent() const noexcept { return np; }

    [[nodiscard]] double plasticity_viscous_denominator() const noexcept { return sp; }

    [[nodiscard]] double damage_viscous_exponent() const noexcept { return nd; }

    [[nodiscard]] double damage_viscous_denominator() const noexcept { return sd; }

protected:
    /// Kinematic hardening numerator
    double gamma = 1.0;
    /// Kinematic hardening denominator
    double C = 1.0;
    /// viscous denominator
    double sp = 1.0;
    /// viscous exponent
    double np = 1.0;
    /// damage viscous denominator
    double sd = 1.0;
    /// damage viscous exponent
    double nd = 1.0;
};
}
