
#pragma once

#include "Tensor.hpp"

namespace neon
{
/**
 * Computes the derivative of the tensor log with respect to \p A
   \f{align*}{
     & \frac{\partial \log \mathbf{A}}{\partial \mathbf{A}}
   \f}
 */
[[nodiscard]] matrix3 log_symmetric_tensor_derivative(matrix2 const& A);
}
