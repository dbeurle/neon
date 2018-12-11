
#include "quadrature/minimum_degree.hpp"

namespace neon
{
auto minimum_degree(int const polynomial_order,
                    int const monomial_order,
                    int const derivative_order) noexcept -> int
{
    return polynomial_order + monomial_order - 2 * derivative_order;
}
}
