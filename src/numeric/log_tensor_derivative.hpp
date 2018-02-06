
#pragma once

#include "numeric/dense_matrix.hpp"

namespace neon
{
/**
 * Computes the derivative of the tensor log with respect to \p A
   \f{align*}{
     & \frac{\partial \log \mathbf{A}}{\partial \mathbf{A}}
   \f}
 * @return Fourth order tensor in Voigt notation
 */
[[nodiscard]] matrix3 log_symmetric_tensor_derivative(matrix2 const& A);
}
