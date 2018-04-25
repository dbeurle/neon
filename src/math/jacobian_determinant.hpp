
#pragma once

#include "numeric/dense_matrix.hpp"
#include <Eigen/Geometry>

/// \file jacobian_determinant.hpp

namespace neon
{
namespace detail
{
template <typename matrix_type>
[[nodiscard]] auto jacobian_determinant(matrix_type const& jacobian)
{
    return jacobian.determinant();
}

[[nodiscard]] inline auto jacobian_determinant(matrix32 const& jacobian)
{
    return jacobian.col(0).cross(jacobian.col(1)).norm();
}
}

/**
 * Compute the Jacobian determinant for a volume or surface mapping.  In three
 * dimensions it performs the following operation:
 * \f{align*}{
 *     j &= \det \begin{bmatrix}
 *             \frac{\partial x}{\partial \xi} & \frac{\partial x}{\partial \eta} & \frac{\partial x}{\partial \zeta} \\
 *             \frac{\partial y}{\partial \xi} & \frac{\partial y}{\partial \eta} & \frac{\partial y}{\partial \zeta} \\
 *             \frac{\partial z}{\partial \xi} & \frac{\partial z}{\partial \eta} & \frac{\partial z}{\partial \zeta}
 *          \end{bmatrix}
 * \f}
 *
 * However the determinant for non-square Jacobian is not defined.  When there
 * is a mapping from \f$ \mathbb{R}^3 \f$ to \f$ \mathbb{R}^2 \f$ the Jacobian
 * `determinant' can be computed by
 * \f{align*}{
 *     j &= || \mathbf{x}_{,\mathbf{\xi}} \times \mathbf{x}_{,\mathbf{\eta}} ||
 * \f}
 * where the Jacobian is given by the non-square matrix
 * \f{align*}{
 *    & \begin{bmatrix}
 *          \frac{\partial x}{\partial \xi} & \frac{\partial x}{\partial \eta} \\
 *          \frac{\partial y}{\partial \xi} & \frac{\partial y}{\partial \eta} \\
 *          \frac{\partial z}{\partial \xi} & \frac{\partial z}{\partial \eta}
 *      \end{bmatrix}
 * \f}
 * This is useful when the surface is described by a two dimensional
 * element but there are only three dimensional coordinates \cite Anton1998.
 */
template <typename matrix_expression>
[[nodiscard]] inline auto jacobian_determinant(matrix_expression const& jacobian)
{
    return detail::jacobian_determinant(jacobian.eval());
}
}
