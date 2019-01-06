
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

    std::vector<IndexType> dof_unique_view(begin(dof_set), end(dof_set));

    for (auto& i : dof_unique_view) i *= dof_size;

    return dof_unique_view;
}

template <typename IndicesType, typename DofOrderType>
void dof_allocator(IndicesType const& node_indices,
                   IndicesType& dof_indices,
                   DofOrderType const dof_order)
{
    // Allocate the degree of freedom indices
    dof_indices.resize(node_indices.rows() * dof_order.size(), node_indices.cols());

    for (typename IndicesType::Index i{0}; i < node_indices.cols(); ++i)
    {
        transform_expand_view(node_indices(Eigen::all, i), dof_indices(Eigen::all, i), dof_order);
    }
}
}
