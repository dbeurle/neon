
#pragma once

#include <type_traits>

namespace neon
{
/// Store the row and column index and return a value of zero.  This is used
/// for computing the sparsity (nonzero) pattern of a sparse matrix using
/// the Eigen linear algebra library
/// \tparam index_type Integer index type
template <class index_type>
class doublet
{
    static_assert(std::is_integral<index_type>::value, "index_type needs to be an integer");

public:
    /// Constructor
    /// \param i Row index
    /// \param j Column index
    doublet(index_type const i, index_type const j) : m_row(i), m_col(j) {}

    /// \return row index of the element
    auto row() const noexcept { return m_row; }

    /// \returns column index of the element
    auto col() const noexcept { return m_col; }

    /// \returns zero
    auto value() const noexcept { return 0.0; }

protected:
    index_type m_row, m_col;
};
}
