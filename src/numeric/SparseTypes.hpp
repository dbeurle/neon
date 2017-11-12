
#pragma once

#include <Eigen/Sparse>

namespace neon
{
using SparseMatrix = Eigen::SparseMatrix<double, Eigen::RowMajor>;

template <class IndexTp>
class Doublet
{
public:
    Doublet(IndexTp const& i, IndexTp const& j) : m_row(i), m_col(j) {}
    /** @return the row index of the element */
    const auto& row() const { return m_row; }
    /** @returns the column index of the element */
    const auto& col() const { return m_col; }
    /** @returns the value of the element */
    auto value() const { return 0.0; }

protected:
    IndexTp m_row, m_col;
};
}
