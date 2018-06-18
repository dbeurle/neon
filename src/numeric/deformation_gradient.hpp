
#pragma once

#include "numeric/dense_matrix.hpp"

/// \file deformation_gradient.hpp

namespace neon::mechanical
{
/// Compute the local deformation gradient
/// \f{align*}{ F_{\xi} &= \bf{x}_\xi \f}
/// \param rhea Shape function gradients at quadrature point
/// \param configuration Configuration of the element (coordinates)
template <typename shape_function_derivatives, typename Configuration>
[[nodiscard]] auto local_deformation_gradient(shape_function_derivatives const& rhea,
                                              Configuration const& configuration)
{
    return configuration * rhea;
}

/// Compute the deformation gradient, F, from the global to local mapping
/// \f{align*}{
/// F &= F_{\xi} \times (F^0_{\xi})^{-1}
/// \f}
/// \param rhea Shape function derivatives at the integration points
/// \param X Reference configuration (spatial coordinates, local nodes)
/// \param x Current configuration (spatial coordinates, local nodes)
[[nodiscard]] inline matrix2 deformation_gradient(matrix const& rhea,
                                                  matrix2x const& X,
                                                  matrix2x const& x)
{
    // Deformation gradient in the reference and current configuration
    return local_deformation_gradient(rhea, x) * local_deformation_gradient(rhea, X).inverse();
}

/// Compute the deformation gradient, F, from the global to local mapping
/// \f{align*}{
/// F &= F_{\xi} \times (F^0_{\xi})^{-1}
/// \f}
/// \param rhea Shape function derivatives at the integration points
/// \param X Reference configuration (spatial coordinates, local nodes)
/// \param x Current configuration (spatial coordinates, local nodes)
[[nodiscard]] inline matrix3 deformation_gradient(matrix const& rhea,
                                                  matrix3x const& X,
                                                  matrix3x const& x)
{
    // Deformation gradient in the reference and current configuration
    return local_deformation_gradient(rhea, x) * local_deformation_gradient(rhea, X).inverse();
}
}
