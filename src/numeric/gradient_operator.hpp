
#pragma once

#include "numeric/dense_matrix.hpp"

#include "constitutive/voigt_dimension.hpp"

namespace neon::fem
{
/**
 * sym_gradient signals the use of
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
template <int spatial_dimension>
inline void sym_gradient(matrix& B, matrix const& L)
{
    auto const nodes_per_element = L.cols();

    if (spatial_dimension == 3)
    {
        for (auto a = 0; a < nodes_per_element; ++a)
        {
            auto const b = a * spatial_dimension;

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
    else if (spatial_dimension == 2)
    {
        for (auto a = 0; a < nodes_per_element; ++a)
        {
            auto const b = a * spatial_dimension;

            B(0, b) = L(0, a);
            B(2, b) = L(1, a);
            B(1, b + 1) = L(1, a);
            B(2, b + 1) = L(0, a);
        }
    }
}

/**
 * sym_gradient signals the use of
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
template <int spatial_dimension>
inline matrix sym_gradient(matrix const& L)
{
    auto const nodes_per_element = L.cols();

    matrix B = matrix::Zero(spatial_to_voigt(spatial_dimension),
                            spatial_dimension * nodes_per_element);

    if (spatial_dimension == 3)
    {
        for (auto a = 0; a < nodes_per_element; ++a)
        {
            auto const b = a * spatial_dimension;

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
    else if (spatial_dimension == 2)
    {
        for (auto a = 0; a < nodes_per_element; ++a)
        {
            auto const b = a * spatial_dimension;

            B(0, b) = L(0, a);
            B(2, b) = L(1, a);
            B(1, b + 1) = L(1, a);
            B(2, b + 1) = L(0, a);
        }
    }
    return B;
}
}
