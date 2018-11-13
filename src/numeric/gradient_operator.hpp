
#pragma once

#include "numeric/dense_matrix.hpp"

#include "constitutive/voigt_dimension.hpp"

/// \file gradient_operator.hpp

namespace neon
{
/** sym_gradient signals the use of
 * \f{align*}{
 * a(w,u) &= (\nabla_S N_a^T) \hspace{0.1cm} (\nabla_S N_a)
 * \f}
 * where in a vector case the symmetric gradient operator
 * takes on a form depending on the number of dimensions.
 *
 * For two dimensions:
 * \f{align*}{
 * \nabla_S &= \begin{bmatrix} N_{a,1} & 0       \\
 *                             0       & N_{a,2} \\
 *                             N_{a,2} & N_{a,1} \\
 * \end{bmatrix}
 * \f}
 * For three dimensions:
 * \f{align*}{
 * \nabla_S &= \begin{bmatrix} N_{a,1} & 0       & 0       \\
 *                             0       & N_{a,2} & 0       \\
 *                             0       &  0      & N_{a,3} \\
 *                             0       & N_{a,3} & N_{a,2} \\
 *                             N_{a,3} & 0       & N_{a,1} \\
 *                             N_{a,2} & N_{a,1} & 0       \\
 * \end{bmatrix}
 * \f}
 * Where \f$ N_{a,i} \f$ represents the partial derivative with respect
 * to the \f$ i \f$th dimension of the shape function in the parent domain
 */
template <int SpatialDimension>
inline void symmetric_gradient(matrix& B, matrix const& L)
{
    auto const nodes_per_element = L.cols();

    if constexpr (SpatialDimension == 3)
    {
        for (auto a = 0; a < nodes_per_element; ++a)
        {
            auto const b = a * SpatialDimension;

            B(0, b) = L(0, a);
            B(4, b) = L(2, a);
            B(5, b) = L(1, a);
            B(1, b + 1) = L(1, a);
            B(3, b + 1) = L(2, a);
            B(5, b + 1) = L(0, a);
            B(2, b + 2) = L(2, a);
            B(3, b + 2) = L(1, a);
            B(4, b + 2) = L(0, a);
        }
    }
    else if (SpatialDimension == 2)
    {
        for (auto a = 0; a < nodes_per_element; ++a)
        {
            auto const b = a * SpatialDimension;

            B(0, b) = L(0, a);
            B(2, b) = L(1, a);
            B(1, b + 1) = L(1, a);
            B(2, b + 1) = L(0, a);
        }
    }
}

template <typename StaticMatrix>
[[nodiscard]] matrix local_gradient(matrix const& dN, StaticMatrix const& Jacobian) noexcept
{
    return (dN * Jacobian.inverse()).transpose();
}
}
