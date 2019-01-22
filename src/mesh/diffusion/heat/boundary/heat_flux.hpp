
#pragma once

#include "interpolations/shape_function.hpp"
#include "mesh/boundary/neumann.hpp"
#include "math/integral_form.hpp"
#include "quadrature/numerical_quadrature.hpp"

namespace neon::diffusion::boundary
{
/// heat_flux is a Neumann type boundary condition
using heat_flux = constant_neumann<fem::integral<surface_interpolation, surface_quadrature, 0>>;
}
