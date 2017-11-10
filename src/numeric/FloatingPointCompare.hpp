
#pragma once

#include <type_traits>

namespace neon
{
template <class T>
typename std::enable_if<std::is_floating_point<T>::value, bool>::type is_approx(T const x,
                                                                                T const y,
                                                                                int const ulp = 2)
{
    // Taken and modified from
    // http://en.cppreference.com/w/cpp/types/numeric_limits/epsilon
    // Since the numeric epsilon is defined at 1.0 then it must be scaled by
    // the worse case (x + y) and accounted for my the ULP (units in the last place).
    return std::abs(x - y) < std::numeric_limits<T>::epsilon() * std::abs(x + y) * ulp
           || std::abs(x - y) < std::numeric_limits<T>::min();
}
}
