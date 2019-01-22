
#pragma once

#include <cmath>
#include <type_traits>

/// \file linear_interpolation.hpp

namespace neon
{
/// relative_distance computes the value in [0, 1] between two end points.
/// It assumes that the value is between the lower and upper values.
/// \param value Inside [lower, upper]
/// \param lower Lower value
/// \param upper Upper value
/// \tparam T Floating point type (float, double, long double)
template <class T>
[[nodiscard]] inline T relative_distance(T const value, T const lower, T const upper) noexcept
{
    static_assert(std::is_floating_point<T>::value, "fraction only accepts floating point arguments");

    return (value - lower) / (upper - lower);
}

/// linear_interpolation performs a linear interpolation using the inbuilt
/// fused-multiply add function \sa std::fma.  This method will be slow unless
/// native optimisations are used but will maintain outstanding accuracy.
/// With the znver1 architecture only two instructions are output
/// \param fraction [0, 1]
/// \param lower Lower value
/// \param upper Upper value
/// \tparam T Floating point type (float, double, long double)
template <class T>
[[nodiscard]] inline T linear_interpolation(T const fraction, T const lower, T const upper)
{
    static_assert(std::is_floating_point<T>::value,
                  "linear_interpolation only accepts floating point arguments");
    return std::fma(fraction, upper, std::fma(-fraction, lower, lower));
}
}
