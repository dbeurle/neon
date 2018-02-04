
#pragma once

#include "DenseMatrix.hpp"

namespace neon
{
/**
 * Performs the tensor dot product on two second order tensors in three dimensions.
 */
template <class MatrixLeft, class MatrixRight>
[[nodiscard]] double double_dot(MatrixLeft const& a, MatrixRight const& b) {
    return (a.array() * b.array()).sum();
}

namespace detail
{
    /** @return the volumetric part of the tensor */
    [[nodiscard]] inline matrix2 volumetric(matrix2 const& a)
    {
        return matrix2::Identity() * a.trace() / 3.0;
    }

    /** @return the volumetric part of the tensor */
    [[nodiscard]] inline matrix3 volumetric(matrix3 const& a)
    {
        return matrix3::Identity() * a.trace() / 3.0;
    }

    /** @return the deviatoric part of the tensor */
    [[nodiscard]] inline matrix2 deviatoric(matrix2 const& a) { return a - volumetric(a); }

    /** @return the deviatoric part of the tensor */
    [[nodiscard]] inline matrix3 deviatoric(matrix3 const& a) { return a - volumetric(a); }
}

/** @return the volumetric part of the tensor */
template <typename MatrixExpression>
[[nodiscard]] auto volumetric(MatrixExpression const& a)
{
    return detail::volumetric(a.eval());
}

/** @return the deviatoric part of the tensor */
template <typename MatrixExpression>
[[nodiscard]] auto deviatoric(MatrixExpression const& a)
{
    return detail::deviatoric(a.eval());
}

/**
 * Evaluates the expression
 * \f{align*}{
   \bar{F} &= \det(F)^{-1/3} F
   \f}
 * @return The unimodular decomposition of a second order tensor
 */
[[nodiscard]] inline matrix3 unimodular(matrix3 const& a)
{
    return std::pow(a.determinant(), -1.0 / 3.0) * a;
}

/** @return The symmetric part of the tensor */
[[nodiscard]] inline matrix2 symmetric(matrix2 const& a) { return 0.5 * (a.transpose() + a); }

/** @return The symmetric part of the tensor */
[[nodiscard]] inline matrix3 symmetric(matrix3 const& a) { return 0.5 * (a.transpose() + a); }

/**
 * I1 returns the coefficient I1, the first stress invariant,
 * which is equal to the trace
 * @return First invariant
 */
[[nodiscard]] inline double I1(matrix3 const& a) { return a.trace(); }

/**
 * I2 returns the coefficient I2, the second stress invariant,
 * which is calculated by:
 * \f{align*}{
 * I_2 &= \frac{1}{2} \left( (tr \tau)^2 - tr(\tau \tau) \right)
 * \f}
 * @return Second invariant
 */
[[nodiscard]] inline double I2(matrix3 const& a)
{
    return 0.5 * (std::pow(a.trace(), 2) - (a * a).trace());
}

/** @return Third invariant, which is the determinant of the tensor */
[[nodiscard]] inline double I3(matrix3 const& a) { return a.determinant(); }

[[nodiscard]] inline matrix identity_expansion(matrix const& H, int const nodal_dofs)
{
    assert(H.rows() == H.cols());
    // Create the geometric part of the tangent stiffness matrix
    matrix K = matrix::Zero(H.rows() * nodal_dofs, H.rows() * nodal_dofs);
    for (auto i = 0; i < H.rows(); ++i)
        for (auto j = 0; j < H.rows(); ++j)
            for (auto k = 0; k < nodal_dofs; ++k)
                K(i * nodal_dofs + k, j * nodal_dofs + k) = H(i, j);
    return K;
}

/*!
 * Handles the representation of common tensors (deviatoric, identity, etc)
 * in Voigt notation suitable for the computation of tensor operations leveraging
 * matrix-vector or matrix-matrix operations.
 * \addtogroup voigt
 * @{
 */
namespace voigt
{
/**
 * Compute the outer product in Voigt notation according to
 * \f$ \mathbb{\mathbf{1} \otimes \mathbf{1}} = \delta_{ij} \delta_{kl} \f$
 */
[[nodiscard]] inline matrix6 I_outer_I()
{
    // clang-format off
    return (matrix6() << 1.0, 1.0, 1.0, 0.0, 0.0, 0.0,
                         1.0, 1.0, 1.0, 0.0, 0.0, 0.0,
                         1.0, 1.0, 1.0, 0.0, 0.0, 0.0,
                         0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                         0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                         0.0, 0.0, 0.0, 0.0, 0.0, 0.0).finished();
    // clang-format on
}

namespace d2
{
/**
 * Compute the outer product in Voigt notation according to
 * \f$ \mathbb{\mathbf{1} \otimes \mathbf{1}} = \delta_{ij} \delta_{kl} \f$
 */
[[nodiscard]] inline matrix3 I_outer_I()
{
    // clang-format off
     return (matrix3() << 1.0, 1.0, 0.0,
                          1.0, 1.0, 0.0,
                          0.0, 0.0, 0.0).finished();
    // clang-format on
}
}

//! Kinematic description of tensor to voigt notation where off diagonal components
//! are multiplied by a factor of two (strain type)
namespace kinematic
{
namespace detail
{
/**
 * Convert second order tensor to Voigt notation according to
 * \f$ \begin{bmatrix} \varepsilon_{11} \\ \varepsilon_{22} \\ 2\varepsilon_{12} \end{bmatrix} =
 * \begin{bmatrix} \varepsilon_{11} & \varepsilon_{12} \\ \varepsilon_{21} &
 * \varepsilon_{22} \end{bmatrix} \f$
 */
[[nodiscard]] inline vector3 to(matrix2 const& a)
{
    // clang-format off
            return (vector3() << a(0, 0),
                                 a(1, 1),
                                 2.0 * a(0, 1)).finished();
    // clang-format on
}

/**
 * Convert second order tensor to Voigt notation according to
 * \f$ \begin{bmatrix} \varepsilon_{11} \\ \varepsilon_{22} \\ \varepsilon_{33} \\ 2\varepsilon_{23}
 * \\ 2\varepsilon_{13} \\
 * 2\varepsilon_{12} \end{bmatrix} = \begin{bmatrix} \varepsilon_{11} & \varepsilon_{12} &
 * \varepsilon_{13} \\ \varepsilon_{21} & \varepsilon_{22} & \varepsilon_{23} \\ \varepsilon_{31} &
 * \varepsilon_{32} & \varepsilon_{33} \end{bmatrix} \f$
 */
[[nodiscard]] inline vector6 to(matrix3 const& a)
{
    // clang-format off
        return (vector6() << a(0, 0),
                             a(1, 1),
                             a(2, 2),
                             2.0 * a(1, 2),
                             2.0 * a(0, 2),
                             2.0 * a(0, 1)).finished();
    // clang-format on
}

/**
 * Convert Voigt notation to second order tensor according to
 * \f$  \begin{bmatrix} \varepsilon_{11} & \varepsilon_{12} \\
 *                      \varepsilon_{21} & \varepsilon_{22}
        \end{bmatrix} = \begin{bmatrix} \varepsilon_{11} \\
                                        \varepsilon_{22} \\
                                        2\varepsilon_{12}
                        \end{bmatrix} \f$
 */
[[nodiscard]] inline matrix2 from(vector3 const& a)
{
    // clang-format off
        return (matrix2() <<     a(0), a(3)/2.0,
                             a(3)/2.0,     a(1)).finished();
    // clang-format on
}

/**
 * Convert Voigt notation to second order tensor according to
 * \f$  \begin{bmatrix} \varepsilon_{11} & \varepsilon_{12} &
 * \varepsilon_{13} \\ \varepsilon_{21} & \varepsilon_{22} & \varepsilon_{23} \\ \varepsilon_{31} &
 * \varepsilon_{32} & \varepsilon_{33} \end{bmatrix} = \begin{bmatrix} \varepsilon_{11} \\
 * \varepsilon_{22} \\ \varepsilon_{33} \\ 2\varepsilon_{23}
 * \\ 2\varepsilon_{13} \\ 2\varepsilon_{12} \end{bmatrix} \f$
 */
[[nodiscard]] inline matrix3 from(vector6 const& a)
{
    // clang-format off
        return (matrix3() <<     a(0), a(5)/2.0, a(4)/2.0,
                             a(5)/2.0,     a(1),     a(3),
                             a(4)/2.0, a(3)/2.0,     a(2)).finished();
    // clang-format on
}
}

template <typename matrix_expression>
[[nodiscard]] inline auto to(matrix_expression const& a)
{
    return detail::to(a.eval());
}

template <typename matrix_expression>
[[nodiscard]] inline auto from(matrix_expression const& a)
{
    return detail::from(a.eval());
}

namespace d2
{
/**
 * Compute the deviatoric tensor in Voigt notation according to
 * \f$ \mathbb{P} = \frac{1}{2}(\delta_{ik} \delta_{jl} + \delta_{il} \delta_{jk}) -
 * \frac{1}{3}\delta_{ij} \delta_{kl} \f$
 */
[[nodiscard]] inline matrix3 deviatoric()
{
    // clang-format off
        return (matrix3() << 2.0 / 3.0, -1.0 / 3.0, 0.0,
                            -1.0 / 3.0,  2.0 / 3.0, 0.0,
                                   0.0,        0.0, 0.5).finished();
    // clang-format on
}
}

/**
 * Compute the deviatoric tensor in Voigt notation according to
 * \f$ \mathbb{P} = \frac{1}{2}(\delta_{ik} \delta_{jl} + \delta_{il} \delta_{jk}) -
 * \frac{1}{3}\delta_{ij} \delta_{kl} \f$
 */
[[nodiscard]] inline matrix6 deviatoric()
{
    // clang-format off
    return (matrix6() << 2.0 / 3.0, -1.0 / 3.0, -1.0 / 3.0, 0.0, 0.0, 0.0,
                        -1.0 / 3.0,  2.0 / 3.0, -1.0 / 3.0, 0.0, 0.0, 0.0,
                        -1.0 / 3.0, -1.0 / 3.0,  2.0 / 3.0, 0.0, 0.0, 0.0,
                               0.0,        0.0,        0.0, 0.5, 0.0, 0.0,
                               0.0,        0.0,        0.0, 0.0, 0.5, 0.0,
                               0.0,        0.0,        0.0, 0.0, 0.0, 0.5).finished();
    // clang-format on
}

/**
 * Compute the fourth order symmetric identity tensor in Voigt notation according to
 * \f$ \mathbb{I} = \frac{1}{2}(\delta_{ik} \delta_{jl} + \delta_{il} \delta_{jk}) \f$
 */
[[nodiscard]] inline matrix6 fourth_order_identity()
{
    // clang-format off
    return (matrix6() << 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                         0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
                         0.0, 0.0, 1.0, 0.0, 0.0, 0.0,
                         0.0, 0.0, 0.0, 0.5, 0.0, 0.0,
                         0.0, 0.0, 0.0, 0.0, 0.5, 0.0,
                         0.0, 0.0, 0.0, 0.0, 0.0, 0.5).finished();
    // clang-format on
}

/**
 * Compute the fourth order symmetric identity tensor in Voigt notation according to
 * \f$ \mathbb{I} = \frac{1}{2}(\delta_{ik} \delta_{jl} + \delta_{il} \delta_{jk}) \f$
 */
[[nodiscard]] inline matrix6 identity()
{
    // clang-format off
    return (matrix6() << 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                         0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
                         0.0, 0.0, 1.0, 0.0, 0.0, 0.0,
                         0.0, 0.0, 0.0, 0.5, 0.0, 0.0,
                         0.0, 0.0, 0.0, 0.0, 0.5, 0.0,
                         0.0, 0.0, 0.0, 0.0, 0.0, 0.5).finished();
    // clang-format on
}

namespace d2
{
/**
 * Compute the fourth order symmetric identity tensor in Voigt notation according to
 * \f$ \mathbb{I} = \frac{1}{2}(\delta_{ik} \delta_{jl} + \delta_{il} \delta_{jk}) \f$
 */
[[nodiscard]] inline matrix3 identity()
{
    // clang-format off
    return (matrix3() << 1.0, 0.0, 0.0,
                         0.0, 1.0, 0.0,
                         0.0, 0.0, 0.5).finished();
    // clang-format on
}
}
}

//! Kinetic description of tensor to voigt notation where off diagonal components
//! are not multiplied by a factor of two
namespace kinetic
{
namespace detail
{
/**
 * Convert second order tensor to Voigt notation according to
 * \f$ \begin{bmatrix} \sigma_{11} \\ \sigma_{22} \\ \sigma_{12} \end{bmatrix} = \begin{bmatrix}
 * \sigma_{11} & \sigma_{12}  \\ \sigma_{21} & \sigma_{22} \end{bmatrix} \f$
 */
[[nodiscard]] inline vector3 to(matrix2 const& a)
{
    return (vector3() << a(0, 0), a(1, 1), a(0, 1)).finished();
}

/**
 * Convert second order tensor to Voigt notation according to
 * \f$ \begin{bmatrix} \sigma_{11} \\ \sigma_{22} \\ \sigma_{33} \\ \sigma_{23} \\ \sigma_{13} \\
 * \sigma_{12} \end{bmatrix} = \begin{bmatrix} \sigma_{11} & \sigma_{12} & \sigma_{13} \\
 * \sigma_{21} & \sigma_{22} & \sigma_{23} \\ \sigma_{31} & \sigma_{32} & \sigma_{33} \end{bmatrix}
 * \f$
 */
[[nodiscard]] inline vector6 to(matrix3 const& a)
{
    return (vector6() << a(0, 0), a(1, 1), a(2, 2), a(1, 2), a(0, 2), a(0, 1)).finished();
}

/**
 * Convert Voigt notation to second order tensor according to
 * \f$ \begin{bmatrix} \sigma_{11} & \sigma_{12} & \sigma_{13} \\
 * \sigma_{21} & \sigma_{22} & \sigma_{23} \\ \sigma_{31} & \sigma_{32} & \sigma_{33} \end{bmatrix}
 * = \begin{bmatrix} \sigma_{11} \\ \sigma_{22} \\ \sigma_{33} \\ \sigma_{23} \\ \sigma_{13} \\
 * \sigma_{12} \end{bmatrix}
 * \f$
 */
[[nodiscard]] inline matrix2 from(vector3 const& a)
{
    // clang-format off
    return (matrix2() << a(0), a(2),
                         a(2), a(1)).finished();
    // clang-format on
}

/**
 * Convert Voigt notation to second order tensor according to
 * \f$ \begin{bmatrix} \sigma_{11} & \sigma_{12} & \sigma_{13} \\
 * \sigma_{21} & \sigma_{22} & \sigma_{23} \\ \sigma_{31} & \sigma_{32} & \sigma_{33} \end{bmatrix}
 * = \begin{bmatrix} \sigma_{11} \\ \sigma_{22} \\ \sigma_{33} \\ \sigma_{23} \\ \sigma_{13} \\
 * \sigma_{12} \end{bmatrix}
 * \f$
 */
[[nodiscard]] inline matrix3 from(vector6 const& a)
{
    // clang-format off
    return (matrix3() << a(0), a(5), a(4),
                         a(5), a(1), a(3),
                         a(4), a(3), a(2)).finished();
    // clang-format on
}
}

template <typename matrix_expression>
[[nodiscard]] inline auto from(matrix_expression const& a)
{
    return detail::from(a.eval());
}

template <typename matrix_expression>
[[nodiscard]] inline auto to(matrix_expression const& a)
{
    return detail::to(a.eval());
}

/**
 * Compute the deviatoric tensor in Voigt notation according to
 * \f$ \mathbb{P} = \frac{1}{2}(\delta_{ik} \delta_{jl} + \delta_{il} \delta_{jk}) -
 * \frac{1}{3}\delta_{ij} \delta_{kl} \f$
 */
[[nodiscard]] inline matrix6 deviatoric()
{
    // clang-format off
    return (matrix6() << 2.0 / 3.0, -1.0 / 3.0, -1.0 / 3.0, 0.0, 0.0, 0.0,
                        -1.0 / 3.0,  2.0 / 3.0, -1.0 / 3.0, 0.0, 0.0, 0.0,
                        -1.0 / 3.0, -1.0 / 3.0,  2.0 / 3.0, 0.0, 0.0, 0.0,
                               0.0,        0.0,        0.0, 1.0, 0.0, 0.0,
                               0.0,        0.0,        0.0, 0.0, 1.0, 0.0,
                               0.0,        0.0,        0.0, 0.0, 0.0, 1.0).finished();
    // clang-format on
}

/**
 * Compute the fourth order symmetric identity tensor in Voigt notation according to
 * \f$ \mathbb{I} = \delta_{ijkl} \f$
 */
[[nodiscard]] inline matrix6 fourth_order_identity() { return matrix6::Identity(); }
}
}
/*! @} End of Doxygen Groups */

/**
 * Computes the outer product of two second order tensors in Voigt notation
 * in two dimensions
    \f{align*}{
        & \mathbf{a} \otimes \mathbf{b} \otimes
    \f}
 * @return fourth order tensor in Voigt notation
 */
[[nodiscard]] inline matrix3 outer_product(vector3 const& a, vector3 const& b)
{
    return a * b.transpose();
}

/**
 * Computes the outer product of two second order tensors in two dimensions
    \f{align*}{
        & \mathbf{a} \otimes \mathbf{b} \otimes
    \f}
 * @return fourth order tensor in Voigt notation
 */
[[nodiscard]] inline matrix3 outer_product(matrix2 const& a, matrix2 const& b)
{
    return voigt::kinetic::to(a) * voigt::kinetic::to(b).transpose();
}

/**
 * Computes the outer product of two second order tensors in three dimensions
    \f{align*}{
        & \mathbf{a} \otimes \mathbf{b} \otimes
    \f}
 * @return fourth order tensor in Voigt notation
 */
[[nodiscard]] inline matrix6 outer_product(matrix3 const& a, matrix3 const& b)
{
    return voigt::kinetic::to(a) * voigt::kinetic::to(b).transpose();
}

/**
 * Computes the outer product of four first order tensors in three dimensions
    \f{align*}{
        & \mathbf{a} \otimes \mathbf{b} \otimes \mathbf{c} \otimes \mathbf{d}
    \f}
*/
[[nodiscard]] inline matrix6 outer_product(vector3 const& a,
                                           vector3 const& b,
                                           vector3 const& c,
                                           vector3 const& d)
{
    return outer_product(outer_product(a, b), outer_product(c, d));
}

/**
 * Computes the outer product of two second order tensors in three dimensions
    \f{align*}{
        & \mathbf{a} \otimes \mathbf{b}
    \f}
*/
[[nodiscard]] inline matrix6 outer_product(matrix3 const& h) { return outer_product(h, h); }

namespace detail
{
[[nodiscard]] inline matrix6 mandel_notation(matrix6 A)
{
    A.block<3, 3>(0, 3) *= std::sqrt(2.0);
    A.block<3, 3>(3, 0) *= std::sqrt(2.0);
    A.block<3, 3>(3, 3) *= 2.0;
    return A;
}

[[nodiscard]] inline matrix3 mandel_notation(matrix3 A)
{
    A.block<1, 2>(2, 0) *= std::sqrt(2.0);
    A.block<2, 1>(0, 2) *= std::sqrt(2.0);
    A(2, 2) *= 2.0;
    return A;
}
}

/**
 * Convert a fourth order tensor in Voigt notation to Mandel notation.  This
 * is useful in computing the double dot product (contraction) between two
 * fourth order tensors using matrix-matrix multiplication.
     \f{align*}{
     \mathbf{c} &= \mathbf{a} : \mathbf{b} \\
      c_{ijkl} &= a_{ijmn} b_{mnkl} \\
      [\mathbf{c}] &= [\mathbf{a}] [\mathbf{b}]
     \f}
 */
template <typename matrix_expression>
[[nodiscard]] inline auto mandel_notation(matrix_expression A)
{
    return detail::mandel_notation(A.eval());
}
}
