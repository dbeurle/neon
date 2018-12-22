
#pragma once

/// @file

#include "quadrature/numerical_quadrature.hpp"

namespace neon
{
template <typename Quadrature>
class shape_function;

using line_interpolation = shape_function<numerical_quadrature<double>>;
using surface_interpolation = shape_function<surface_quadrature>;
using volume_interpolation = shape_function<volume_quadrature>;
}
