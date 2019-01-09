
#pragma once

#include "quadrature/numerical_quadrature.hpp"

namespace neon
{
/// \cite Hughes2012
class gauss_legendre_line : public line_quadrature
{
public:
    ///
    explicit gauss_legendre_line(int const minimum_degree) noexcept(false);
};
}
