
#pragma once

#include "numeric/index_types.hpp"

#include <algorithm>
#include <set>
#include <vector>

/// \file unique_dof_allocator.hpp

namespace neon
{
/// Allocate a vector of indices containing the non-repeated DoFs from the
/// collection of boundary meshes provided as a function argument.
/// \tparam dof_size Number of DoFs
/// \param boundary_meshes A vector a boundary meshes to extract DoFs from
template <int dof_size, typename boundary_mesh_type, class index_type = std::int32_t>
[[nodiscard]] auto unique_dof_allocator(boundary_mesh_type const& boundary_meshes)
{
    std::set<index_type> dof_set;

    for (auto const& boundary_mesh : boundary_meshes)
    {
        auto const unique_view = boundary_mesh.unique_connectivity();

        std::copy(begin(unique_view), end(unique_view), std::inserter(dof_set, end(dof_set)));
    }

    std::vector<index_type> dof_unique_view(begin(dof_set), end(dof_set));

    for (auto& i : dof_unique_view) i *= dof_size;

    return dof_unique_view;
}
}
