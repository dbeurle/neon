
#pragma once

#include "numeric/DenseMatrix.hpp"

namespace neon
{
/**
 * Perform the spectral decomposition of the provided matrix.  The returned
 * results are:
 * - boolean flag (true if the values are unique)
 * - a pair containing eigen values
 * - a pair containing eigen projections
 */
[[nodiscard]] std::tuple<bool, std::pair<double, double>, std::pair<matrix2, matrix2>> spectral_decomposition(
    matrix2 const& A);
}
