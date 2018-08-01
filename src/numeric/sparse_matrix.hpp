
#pragma once

#define EIGEN_SPARSEMATRIX_PLUGIN "numeric/thread_safe_plugin.hpp"

#include <Eigen/Sparse>

#include <complex>

namespace neon
{
/// Type alias for row major real-valued sparse matrix
using sparse_matrix = Eigen::SparseMatrix<double, Eigen::RowMajor>;
/// Type alias for row major complex-valued sparse matrix
using complex_sparse_matrix = Eigen::SparseMatrix<std::complex<double>, Eigen::RowMajor>;
}
