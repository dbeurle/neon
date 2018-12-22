
#pragma once

/// @file

#include "mesh/boundary/neumann.hpp"
#include "interpolations/shape_function.hpp"

namespace neon::mechanics::plane
{
/// traction is a non-follower load that has a surface interpolation and
/// computes the element external load vector contribution
/// \sa nonfollower_load_boundary
using traction = surface_load<line_interpolation>;
}
