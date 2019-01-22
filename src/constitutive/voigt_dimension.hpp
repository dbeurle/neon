
#pragma once

/// @file

namespace neon
{
/// Convert a spatial dimension to the size of the symmetric Voigt notation
[[deprecated]] constexpr auto spatial_to_voigt(int const spatial_dimension)
{
    return spatial_dimension == 1 ? 1 : spatial_dimension == 2 ? 3 : 6;
}
}
