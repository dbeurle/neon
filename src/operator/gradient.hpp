
#pragma once

namespace neon
{
namespace eulerian
{
/** gradient is defined as
 * \f{align*}{
 * a(w,u) &= (\nabla_S N_a^T) \hspace{0.1cm} (\nabla_S N_a)
 * \f}
 * where in a vector case the symmetric gradient operator
 * takes on a form depending on the number of dimensions.
 * \f{align*}{
 * \nabla_S &= \begin{bmatrix} N_{a,1} & 0       & 0       \\
 *                             0       & N_{a,2} & 0       \\
 *                             0       &  0      & N_{a,3} \\
 *                             0       & N_{a,3} & N_{a,2} \\
 *                             N_{a,3} & 0       & N_{a,1} \\
 *                             N_{a,2} & N_{a,1} & 0       \\
 * \end{bmatrix}
 * \f}
 * where \f$ N_{a,i} \f$ represents the partial derivative with respect
 * to the \f$ i \f$th dimension of the shape function in the parent domain
 */
auto gradient = [](matrixdx<6>& gradient_operator, matrixdx<3> const& local_gradient) {
    auto const nodes_per_element = local_gradient.cols();

    for (auto a = 0; a < nodes_per_element; ++a)
    {
        auto const b = a * 3;

        gradient_operator(0, b) = local_gradient(0, a);
        gradient_operator(4, b + 2) = local_gradient(0, a);
        gradient_operator(5, b + 1) = local_gradient(0, a);
        gradient_operator(5, b) = local_gradient(1, a);
        gradient_operator(1, b + 1) = local_gradient(1, a);
        gradient_operator(3, b + 2) = local_gradient(1, a);
        gradient_operator(2, b + 2) = local_gradient(2, a);
        gradient_operator(3, b + 1) = local_gradient(2, a);
        gradient_operator(4, b) = local_gradient(2, a);
    }
};

namespace plane
{
/** gradient is defined as
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
 * where \f$ N_{a,i} \f$ represents the partial derivative with respect
 * to the \f$ i \f$th dimension of the shape function in the parent domain
 */
auto gradient = [](matrixdx<3>& gradient_operator, matrixdx<2> const& local_gradient) {
    for (auto a = 0; a < nodes_per_element; ++a)
    {
        auto const b = a * 2;

        gradient_operator(0, b) = local_gradient(0, a);
        gradient_operator(2, b + 1) = local_gradient(0, a);
        gradient_operator(2, b) = local_gradient(1, a);
        gradient_operator(1, b + 1) = local_gradient(1, a);
    }
};
}

namespace axisymmetric
{
/// Torsionless symmetric gradient operator for axisymmetric problems
auto gradient = [](matrix4d& gradient_operator,
                   matrix2d& local_gradient,
                   vector const& local_function,
                   double const radius) {
    for (auto a = 0; a < local_gradient.cols(); ++a)
    {
        auto const b = a * 2;

        gradient_operator(0, b) = gradient_operator(2, b + 1) = local_gradient(0, a);
        gradient_operator(2, b) = gradient_operator(1, b + 1) = local_gradient(1, a);
        gradient_operator(3, b) = local_function(a) / radius;
    }
};
}

namespace lagrangian
{
}

}
