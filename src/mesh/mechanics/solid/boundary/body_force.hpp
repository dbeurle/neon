
#pragma once

#include "math/integral_form.hpp"
#include "interpolations/shape_function.hpp"
#include "quadrature/numerical_quadrature.hpp"

#include "mesh/boundary/neumann.hpp"

namespace neon::mechanics::solid
{
/// body_force is a non-follower load that has a volume interpolation and
/// computes the element external load vector contribution
/// \sa NonFollowerLoad
using body_force = constant_neumann<fem::integral<volume_interpolation, volume_quadrature, 0>>;
}
