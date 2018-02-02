
#include "spectral_decomposition.hpp"

#include "numeric/FloatingPointCompare.hpp"

namespace neon
{
std::tuple<bool, std::pair<double, double>, std::pair<matrix2, matrix2>> spectral_decomposition(
    matrix2 const& A)
{
    // First invariant
    auto const I1 = A.trace();
    // Second invariant
    auto const I2 = A.determinant();

    // First eigenvalue
    auto const x1 = (I1 + std::sqrt(std::pow(I1, 2) - 4.0 * I2)) / 2.0;
    // Second eigenvalue
    auto const x2 = (I1 - std::sqrt(std::pow(I1, 2) - 4.0 * I2)) / 2.0;

    // Eigenprojections
    return is_approx(x1, x2)
               ? std::tuple(false,
                            std::pair(x1, x2),
                            std::pair(matrix2::Identity().eval(), matrix2::Identity().eval()))
               : std::tuple(true,
                            std::pair(x1, x2),
                            std::pair(1.0 / (2.0 * x1 - I1)
                                          * (A + (x1 - I1) * matrix2::Identity()).eval(),
                                      1.0 / (2.0 * x2 - I1)
                                          * (A + (x2 - I1) * matrix2::Identity()).eval()));
}
}
