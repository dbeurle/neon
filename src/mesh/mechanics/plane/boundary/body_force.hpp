
#pragma once

#include "mesh/boundary/neumann.hpp"
#include "interpolations/shape_function.hpp"

namespace neon::mechanics::plane
{
/// body_force is a non-follower load that has a volume interpolation and
/// computes the element external load vector contribution
/// \sa nonfollower_load_boundary
using body_force = volume_load<surface_interpolation>;
}
