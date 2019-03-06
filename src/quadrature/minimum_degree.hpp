
#pragma once

#include <cassert>

// \file minumum_degree.hpp

namespace neon
{
/// Compute the minimum degree for numerical quadrature to 'exactly' integrate
/// a given bilinear or linear form
/// \param polynomial_order The highest polynomial order in the expression
/// \param monomial_order The highest monomial order in the expression
/// \param derivative_order Order of derivative in weak form (1 for elasticity, 2 for EB beam theory)
constexpr auto minimum_degree(int const polynomial_order,
                              int const monomial_order,
                              int const derivative_order) noexcept -> int
{
    assert(polynomial_order >= 0);
    assert(monomial_order >= 0);
    assert(derivative_order >= 0);

    return polynomial_order + monomial_order - 2 * derivative_order;
}
}
