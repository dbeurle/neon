
#pragma once

#include "DenseTypes.hpp"

namespace neon
{
namespace fem
{
/**
 * GalerkinSymGradient signals the use of
 * \f{align*}{
 * a(w,u) &= (\nabla_S N_a^T) \hspace{0.1cm} (\nabla_S N_a)
 * \f}
 * Where in a vector case the symmetric gradient operator
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
template <int Dimension>
inline Matrix SymGradient(Matrix const& L)
{
    auto nodesPerElement = L.cols();

    Matrix B = Matrix::Zero(Dimension == 3 ? 6 : 3, Dimension * nodesPerElement);

    // These if statements should be optimized away
    if (Dimension == 3)
    {
        for (auto a = 0; a < nodesPerElement; ++a)
        {
            auto const b = a * Dimension;

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
    else if (Dimension == 2)
    {
        for (int a = 0; a < nodesPerElement; ++a)
        {
            const int b = a * Dimension;
            B(0, b) = L(0, a);
            B(2, b) = L(1, a);
            B(1, b + 1) = L(1, a);
            B(2, b + 1) = L(0, a);
        }
    }
    return B;
}
}
}
