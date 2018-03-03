
#pragma once

#include "mesh/generic/neumann.hpp"

#include "interpolations/shape_function.hpp"

namespace neon::mechanical::solid
{
/**
 * Traction is a non-follower load that has a surface interpolation and
 * computes the element external load vector contribution to the system of
 * equations
 * \sa NonFollowerLoad
 */
using traction = surface_load<surface_interpolation>;
}
