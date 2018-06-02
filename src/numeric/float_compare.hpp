
#pragma once

#include <cmath>
#include <limits>
#include <type_traits>

/// \file float_compare.hpp

/// neon namespace
namespace neon
{
/// Perform a floating point comparison using a specified number of units
/// in the last place
/// \tparam Floating point type
/// \param x Value 1
/// \param y Value 2
/// \param ulp Number of units in last place
template <class T>
std::enable_if_t<std::is_floating_point<T>::value, bool> is_approx(T const x,
                                                                   T const y,
                                                                   int const ulp = 2) noexcept
{
    // Taken and modified from
    // http://en.cppreference.com/w/cpp/types/numeric_limits/epsilon
    // Since the numeric epsilon is defined at 1.0 then it must be scaled by
    // the worse case (x + y) and accounted for by the ULP (units in the last place).
    return std::abs(x - y) < std::numeric_limits<T>::epsilon() * std::abs(x + y) * ulp
           || std::abs(x - y) < std::numeric_limits<T>::min();
}
}
