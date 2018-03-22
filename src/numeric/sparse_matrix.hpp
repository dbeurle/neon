
#pragma once

#define EIGEN_SPARSEMATRIX_PLUGIN "numeric/thread_safe_plugin.hpp"

#include <Eigen/Sparse>

#include <type_traits>

namespace neon
{
/// Type alias for the sparse matrix type
using sparse_matrix = Eigen::SparseMatrix<double, Eigen::RowMajor>;

template <class index_type>
class Doublet
{
    static_assert(std::is_integral<index_type>::value, "index_type needs to be an integer");

public:
    Doublet(index_type const i, index_type const j) : m_row(i), m_col(j) {}
    /// \return the row index of the element
    auto row() const noexcept { return m_row; }
    /// \returns the column index of the element
    auto col() const noexcept { return m_col; }
    /// \returns the value of the element
    auto value() const noexcept { return 0.0; }

protected:
    index_type m_row, m_col;
};
}
