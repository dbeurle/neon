
#pragma once

#include "interpolations/ShapeFunction.hpp"
#include "mesh/generic/Neumann.hpp"

namespace neon::mechanical::solid
{
/**
 * BodyForce is a non-follower load that has a volume interpolation and
 * computes the element external load vector contribution to the system of
 * equations
 * \sa NonFollowerLoad
 */
using body_force = VolumeLoad<VolumeInterpolation>;
}
