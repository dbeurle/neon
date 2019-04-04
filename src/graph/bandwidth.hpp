
#pragma once

#include <algorithm>

namespace neon
{
template <class SparseMatrix>
auto compute_bandwidth(SparseMatrix const& A) noexcept
{
    using index_type = typename SparseMatrix::StorageIndex;

    index_type bandwidth = 0;

    for (index_type row_index = 0; row_index < A.innerSize() - 1; ++row_index)
    {
        bandwidth = std::max(bandwidth,
                             A.outerIndexPtr()[A.innerIndexPtr()[row_index + 1] - 1] - row_index);
    }
    return bandwidth;
}
}
