
#pragma once

#define EIGEN_SPARSEMATRIX_PLUGIN "numeric/thread_safe_plugin.hpp"

#include <Eigen/Sparse>

#include "numeric/doublet.hpp"

namespace neon
{
/// Type alias for the row major sparse matrix type
using sparse_matrix = Eigen::SparseMatrix<double, Eigen::RowMajor>;
}
