
#pragma once

#include "DenseTypes.hpp"
#include <Eigen/SparseCore>

namespace neon
{
using SparseMatrix = Eigen::SparseMatrix<double>;

class Doublet
{
public:
    Doublet(const int& i, const int& j) : m_row(i), m_col(j) {}
    /** @return the row index of the element */
    const auto& row() const { return m_row; }
    /** @returns the column index of the element */
    const auto& col() const { return m_col; }
    /** @returns the value of the element */
    auto value() const { return 0.0; }

protected:
    int m_row, m_col;
};
}
