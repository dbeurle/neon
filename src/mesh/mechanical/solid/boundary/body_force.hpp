
#pragma once

#include "interpolations/shape_function.hpp"
#include "mesh/generic/neumann.hpp"

namespace neon::mechanical::solid
{
/**
 * BodyForce is a non-follower load that has a volume interpolation and
 * computes the element external load vector contribution to the system of
 * equations
 * \sa NonFollowerLoad
 */
using body_force = volume_load<volume_interpolation>;
}
