
#pragma once

#include "numeric/sparse_matrix.hpp"

#include <vector>

namespace neon
{
class adjacency_graph
{
public:
    using index_type = sparse_matrix::StorageIndex;

public:
    explicit adjacency_graph(sparse_matrix const& A) { this->update(A); }

    void update(sparse_matrix const& A)
    {
        m_adjacencies.resize(A.rows());

        for (index_type index = 0; index < A.rows(); ++index)
        {
            auto& adjacency_list = m_adjacencies[index];

            for (sparse_matrix::InnerIterator iterator(A, index); iterator; ++iterator)
            {
                adjacency_list.emplace_back(A.IsRowMajor ? iterator.col() : iterator.row());
            }

            auto const location = std::find(begin(adjacency_list), end(adjacency_list), index);

            if (location != end(adjacency_list))
            {
                adjacency_list.erase(location);
            }
        }
    }

    auto size() const noexcept -> std::size_t { return m_adjacencies.size(); }

    auto children(index_type const node) const noexcept -> std::vector<index_type> const&
    {
        return m_adjacencies[node];
    }

    auto children(index_type const node) noexcept -> std::vector<index_type>&
    {
        return m_adjacencies[node];
    }

    auto degree(index_type const node) const noexcept -> index_type
    {
        return m_adjacencies[node].size();
    }

    void minimum_degree_sort() noexcept(false)
    {
        for (auto& adjacency_list : m_adjacencies)
        {
            std::sort(begin(adjacency_list),
                      end(adjacency_list),
                      [&](index_type const& left, index_type const& right) {
                          return degree(left) < degree(right);
                      });
        }
    }

protected:
    std::vector<std::vector<index_type>> m_adjacencies;
};
}
