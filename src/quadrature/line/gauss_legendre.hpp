
#pragma once

#include "quadrature/numerical_quadrature.hpp"

namespace neon::quadrature::line
{
/// \cite Hughes2012
class gauss_legendre : public line_quadrature
{
public:
    ///
    explicit gauss_legendre(int const minimum_degree) noexcept(false);
};
}
