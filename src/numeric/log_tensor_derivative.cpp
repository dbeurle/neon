
#include "numeric/log_tensor_derivative.hpp"

#include "numeric/spectral_decomposition.hpp"
#include "numeric/tensor_operations.hpp"

namespace neon
{
matrix3 log_symmetric_tensor_derivative(matrix2 const& A)
{
    auto const[is_unique, eigenvalues, eigenprojections] = spectral_decomposition(A);

    auto const[x1, x2] = eigenvalues;

    auto const Isym = voigt::kinematic::d2::identity();

    if (is_unique)
    {
        auto const[E1, E2] = eigenprojections;

        auto const E1_outer_E1 = outer_product(voigt::kinematic::to(E1), voigt::kinematic::to(E1));
        auto const E2_outer_E2 = outer_product(voigt::kinematic::to(E2), voigt::kinematic::to(E2));

        // Derivative of log(x) is 1/x
        auto const y1 = 1.0 / x1;
        auto const y2 = 1.0 / x2;

        return (y2 - y1) / (x2 - x1) * (Isym - E1_outer_E1 - E2_outer_E2) + E1_outer_E1 / x1
               + E2_outer_E2 / x2;
    }
    return Isym / x1;
}
}
