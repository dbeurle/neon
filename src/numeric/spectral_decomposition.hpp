
#pragma once

/// @file

#include "numeric/dense_matrix.hpp"

/// \brief Spectral decomposition routines for small matrices
namespace neon
{
/// Perform the spectral decomposition of the provided matrix.  The returned
/// results are:
/// - boolean flag (true if the values are unique)
/// - a pair containing eigenvalues
/// - a pair containing eigenprojections
[[nodiscard]] std::tuple<bool, std::pair<double, double>, std::pair<matrix2, matrix2>> spectral_decomposition(
    matrix2 const& A);
}
