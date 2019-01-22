
#pragma once

#include "quadrature/numerical_quadrature.hpp"

namespace neon::quadrature::triangle
{
class cowper : public surface_quadrature
{
public:
    explicit cowper(int const minimum_degree) noexcept(false);
};
}
