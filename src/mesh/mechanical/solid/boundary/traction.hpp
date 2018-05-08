
#pragma once

#include "mesh/generic/neumann.hpp"

#include "interpolations/shape_function.hpp"

namespace neon::mechanical::solid
{
/// traction is a non-follower load that has a surface interpolation and
/// computes the element external load vector contribution
/// \sa NonFollowerLoad
using traction = surface_load<surface_interpolation>;
}
