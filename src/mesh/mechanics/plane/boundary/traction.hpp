
#pragma once

#include "mesh/boundary/neumann.hpp"

#include "math/integral_form.hpp"
#include "quadrature/numerical_quadrature.hpp"
#include "interpolations/shape_function.hpp"

namespace neon::mechanics::plane
{
/// traction is a non-follower load that has a surface interpolation and
/// computes the element external load vector contribution
/// \sa nonfollower_load_boundary
using traction = constant_neumann<fem::integral<line_interpolation, line_quadrature, 0>>;
}
