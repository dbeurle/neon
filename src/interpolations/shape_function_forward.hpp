
#pragma once

namespace neon
{
template <typename... Spaces>
class shape_function;

using line_interpolation = shape_function<double>;
using surface_interpolation = shape_function<double, double>;
using volume_interpolation = shape_function<double, double, double>;
}
