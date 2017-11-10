
#pragma once

#include <type_traits>

namespace neon
{
/** Convert a spatial dimension to the size of the symmetric Voigt notation */
constexpr auto spatial_to_voigt(typename std::integral_constant<int, 1>) { return 1; }
constexpr auto spatial_to_voigt(typename std::integral_constant<int, 2>) { return 3; }
constexpr auto spatial_to_voigt(typename std::integral_constant<int, 3>) { return 6; }
}
