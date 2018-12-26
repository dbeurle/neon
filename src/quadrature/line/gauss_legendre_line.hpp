
#pragma once

#include "quadrature/numerical_quadrature.hpp"

namespace neon
{
class gauss_legendre_line : public numerical_quadrature<double>
{
public:
    ///
    explicit gauss_legendre_line(int const minimum_degree) noexcept(false);
};
}
