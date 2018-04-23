
#pragma once

#include "numeric/index_types.hpp"

#include <algorithm>
#include <set>
#include <vector>

namespace neon
{
template <std::int8_t dof_size, typename boundary_mesh_type>
[[nodiscard]] std::vector<std::int32_t> mesh_dof_filter(boundary_mesh_type const& boundary_meshes) {
    std::set<std::int32_t> dof_set;

    for (auto const& boundary_mesh : boundary_meshes)
    {
        auto const unique_view = boundary_mesh.unique_connectivity();

        std::copy(std::begin(unique_view),
                  std::end(unique_view),
                  std::inserter(dof_set, std::end(dof_set)));
    }

    std::vector<std::int32_t> dof_unique_view(std::begin(dof_set), std::end(dof_set));

    for (auto& i : dof_unique_view) i *= dof_size;

    return dof_unique_view;
}
}
