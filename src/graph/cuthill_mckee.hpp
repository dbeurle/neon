
#pragma once

#include "graph/adjacency_graph.hpp"

#include <queue>
#include <vector>

namespace neon
{
/// Find the first lowest degree node in the graph
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
void update_queue(std::queue<T>& queue,
                  std::vector<T> const& children,
                  std::vector<bool> const& is_reordered) noexcept(false)
{
    for (auto const& child : children)
    {
        if (is_reordered[child]) continue;

        queue.push(child);
    }
}

class reverse_cuthill_mcgee
{
public:
    using index_type = sparse_matrix::StorageIndex;

public:
    explicit reverse_cuthill_mcgee(sparse_matrix const& A) : graph(A) {}

    void update(sparse_matrix const& A) noexcept(false) { graph.update(A); }

    void compute() noexcept(false)
    {
        graph.minimum_degree_sort();

        m_permutation.clear();
        m_permutation.reserve(graph.size());

        std::vector<bool> is_reordered(graph.size(), false);

        std::queue<index_type> queue;

        while (m_permutation.size() != graph.size())
        {
            auto const peripheral = find_lowest_degree(is_reordered, graph);

            m_permutation.emplace_back(peripheral);

            is_reordered[peripheral] = true;

            update_queue(queue, graph.children(peripheral), is_reordered);

            while (!queue.empty())
            {
                if (!is_reordered[queue.front()])
                {
                    is_reordered[queue.front()] = true;
                    m_permutation.emplace_back(queue.front());

                    update_queue(queue, graph.children(queue.front()), is_reordered);
                }
                queue.pop();
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
