
#pragma once

/// @file

#include "quadrature/numerical_quadrature.hpp"

namespace neon::quadrature::quadrilateral
{
class gauss_legendre : public surface_quadrature
{
public:
    explicit gauss_legendre(int const minimum_degree) noexcept(false);
};
}
