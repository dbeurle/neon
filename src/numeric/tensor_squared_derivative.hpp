
#pragma once

#include "numeric/tensor_operations.hpp"

namespace neon
{
/**
 * Computes the derivative of the tensor squared with respect to the
 * same tensor
   \f{align*}{
    \mathbf{D}(\mathbf{X}) &= \frac{\partial \mathbf{X}^2}{\partial \mathbf{X}}
   \f}
 */
inline matrix6 tensor_squared_derivative(matrix3 const& X)
{
    matrix6 dx2_dx;
    // clang-format off
    dx2_dx << 2*X(0, 0),       0.0,         0.0,               0.0,           X(2, 0), X(1, 0),
                    0.0, 2*X(1, 1),         0.0,           X(2, 1),               0.0, X(1, 0),
                    0.0,       0.0, 2.0*X(2, 2),           X(2, 1),           X(2, 0), 0.0,
                    0.0,   X(2, 1),     X(2, 1), X(1, 1) + X(2, 2),           X(1, 0), 0.0,
                X(2, 0),         0,     X(2, 0),           X(1, 0), X(0, 0) + X(2, 2), X(1, 2),
                X(1, 0),   X(1, 0),         0.0,               0.0,           X(1, 2), X(0, 0) + X(1, 1);
    // clang-format on
    return dx2_dx;
}
}
