
#pragma once

#include "interpolations/ShapeFunction.hpp"
#include "mesh/generic/neumann.hpp"

namespace neon::mechanical::solid
{
/**
 * BodyForce is a non-follower load that has a volume interpolation and
 * computes the element external load vector contribution to the system of
 * equations
 * \sa NonFollowerLoad
 */
using body_force = boundary::volume_load<VolumeInterpolation>;
}
