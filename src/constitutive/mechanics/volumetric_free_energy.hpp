
#pragma once

/// \file volumetric_free_energy.hpp

namespace neon
{
/**
 * Compute the first derivative of the volumetric free energy function
 * \f{align*}{
     U' &= \frac{\partial U}{\partial J} = \frac{K}{2}\left(J - \frac{1}{J}\right)
   \f}
 * where
 * \f$ U = \frac{K}{4}(J^2 - 1) - \frac{K}{2}\ln{J} \f$
 */
inline double volumetric_free_energy_dJ(double const J, double const bulk_modulus) noexcept
{
    return bulk_modulus / 2.0 * (J - 1.0 / J);
}

/**
 * Compute the second derivative of the volumetric free energy function
 * \f{align*}{
     U^{''} &= \frac{\partial^2 U}{\partial J^2} = \frac{K}{2} \left(1 + \frac{1}{J^2}\right) \f}
 * \sa volumetric_free_energy_dJ
 */
inline double volumetric_free_energy_second_d2J(double const J, double const bulk_modulus) noexcept
{
    return bulk_modulus / 2.0 * (1.0 + 1.0 / std::pow(J, 2));
}
}
