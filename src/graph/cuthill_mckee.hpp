
#pragma once

#include "graph/adjacency_graph.hpp"

#include <stack>
#include <vector>

namespace neon
{
/// Find the first lowest degree node in the graph
/// \tparam
inline auto find_lowest_degree(std::vector<bool> const& is_reordered,
                               adjacency_graph const& graph) noexcept -> adjacency_graph::index_type
{
    using index_type = adjacency_graph::index_type;

    index_type count = graph.size(), node = graph.size();

    for (index_type row = 0, rows = graph.size(); row < rows; ++row)
    {
        auto const row_adjacency = graph.degree(row);

        if (row_adjacency < count && !is_reordered[row])
        {
            count = row_adjacency;
            node = row;
        }
    }
    return node;
}

template <class T>
void update_queue(std::stack<T>& queue,
                  std::vector<T> const& children,
                  std::vector<bool> const& is_reordered) noexcept(false)
{
    for (auto const& child : children)
    {
        if (!is_reordered[child])
        {
            queue.push(child);
        }
    }
}

class reverse_cuthill_mcgee
{
public:
    using index_type = sparse_matrix::StorageIndex;

public:
    reverse_cuthill_mcgee(sparse_matrix const& A) : graph(A) {}

    void update(sparse_matrix const& A) noexcept(false) { graph.update(A); }

    void compute() noexcept(false)
    {
        graph.minimum_degree_sort();

        m_permutation.clear();
        m_permutation.reserve(graph.size());

        std::vector<bool> is_reordered(graph.size(), false);

        std::stack<index_type> queue;

        while (m_permutation.size() != graph.size())
        {
            auto const parent = find_lowest_degree(is_reordered, graph);

            m_permutation.emplace_back(parent);

            is_reordered[parent] = true;

            update_queue(queue, graph.children(parent), is_reordered);

            while (!queue.empty())
            {
                if (!is_reordered[queue.top()])
                {
                    is_reordered[m_permutation.emplace_back(queue.top())] = true;

                    auto const& children = graph.children(queue.top());

                    queue.pop();

                    update_queue(queue, children, is_reordered);
                }
                else
                {
                    queue.pop();
                }
            }
        }
        // Reverse Cuthill-McKee algorithm
        std::reverse(begin(m_permutation), end(m_permutation));
    }

    auto permutation() const noexcept -> std::vector<index_type> const& { return m_permutation; }

protected:
    std::vector<index_type> m_permutation;
    adjacency_graph graph;
};

}
