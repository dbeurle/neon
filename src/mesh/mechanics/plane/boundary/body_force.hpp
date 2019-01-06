
#pragma once

#include "mesh/boundary/neumann.hpp"

#include "quadrature/numerical_quadrature.hpp"
#include "interpolations/shape_function.hpp"
#include "math/integral_form.hpp"

namespace neon::mechanics::plane
{
/// body_force is a non-follower load that has a volume interpolation and
/// computes the element external load vector contribution
/// \sa nonfollower_load_boundary
using body_force = constant_neumann<fem::integral<surface_interpolation, surface_quadrature, 0>>;
}
