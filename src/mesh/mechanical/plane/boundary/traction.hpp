
#pragma once

#include "mesh/generic/neumann.hpp"
#include "interpolations/shape_function.hpp"

namespace neon::mechanical::plane
{
/// traction is a non-follower load that has a surface interpolation and
/// computes the element external load vector contribution
/// \sa nonfollower_load_boundary
using traction = surface_load<line_interpolation>;
}
