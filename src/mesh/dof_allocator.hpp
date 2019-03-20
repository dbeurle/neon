
#pragma once

/// @file

#include "numeric/index_types.hpp"
#include "math/view.hpp"

#include <algorithm>
#include <set>
#include <vector>

namespace neon
{
/// Allocate a vector of indices containing the non-repeated DoFs from the
/// collection of boundary meshes provided as a function argument.
/// \tparam dof_size Number of DoFs
/// \param boundary_meshes A vector a boundary meshes to extract DoFs from
template <int dof_size, typename BoundaryMeshType, class IndexType = std::int32_t>
[[nodiscard]] auto unique_dof_allocator(BoundaryMeshType const& boundary_meshes)
{
    std::set<IndexType> dof_set;

    for (auto const& boundary_mesh : boundary_meshes)
    {
        auto const unique_view = boundary_mesh.unique_node_indices();
        dof_set.insert(begin(unique_view), end(unique_view));
    }

    std::vector<IndexType> unique_indices(dof_set.size());

    std::transform(begin(dof_set), end(dof_set), begin(unique_indices), [=](auto const index) {
        return index * dof_size;
    });

    return unique_indices;
}

/// Allocate the degree of freedom indices
template <typename IndicesType>
void dof_allocator(IndicesType const& node_indices,
                   IndicesType& dof_indices,
                   std::int64_t const dof_order)
{
    dof_indices.resize(node_indices.rows() * dof_order, node_indices.cols());
    transform_expand_n(node_indices.data(), node_indices.size(), dof_indices.data(), dof_order);
}
}
