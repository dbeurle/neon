
#pragma once

/// @file

#include "quadrature/numerical_quadrature.hpp"

namespace neon::quadrature::hexahedron
{
class gauss_legendre : public volume_quadrature
{
public:
    explicit gauss_legendre(int const minimum_degree) noexcept(false);
};
}
