
#pragma once

/// @file

#include <type_traits>

#include <boost/math/special_functions/relative_difference.hpp>

/// neon namespace
namespace neon
{
/// Perform a floating point comparison using boost for large numbers and an absolute error for smaller ones

/// \tparam Floating point type
/// \param x Value 1
/// \param y Value 2
/// \param absolute_tolerance
/// \param relative_tolerance
template <class T>
std::enable_if_t<std::is_floating_point<T>::value, bool> is_approx(T const x,
                                                                   T const y,
                                                                   T const absolute_tolerance = 0.0,
                                                                   T const relative_tolerance = 1e-10) noexcept
{
    // References:
    // https://www.boost.org/doc/libs/1_67_0/libs/math/doc/html/math_toolkit/float_comparison.html
    // https://www.python.org/dev/peps/pep-0485/
    return boost::math::relative_difference(x, y) <= relative_tolerance
           || std::abs(x - y) <= absolute_tolerance;
}
}
