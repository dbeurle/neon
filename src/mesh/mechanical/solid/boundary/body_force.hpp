
#pragma once

#include "interpolations/shape_function.hpp"
#include "mesh/boundary/neumann.hpp"

namespace neon::mechanical::solid
{
/// body_force is a non-follower load that has a volume interpolation and
/// computes the element external load vector contribution
/// \sa NonFollowerLoad
using body_force = volume_load<volume_interpolation>;
}
