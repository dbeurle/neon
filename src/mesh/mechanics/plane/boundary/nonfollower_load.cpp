
#include "nonfollower_load.hpp"

#include "interpolations/interpolation_factory.hpp"

#include "io/json.hpp"

namespace neon::mechanics::plane
{
nonfollower_load_boundary::nonfollower_load_boundary(
    std::shared_ptr<material_coordinates>& coordinates,
    std::vector<basic_submesh> const& submeshes,
    json const& simulation_data,
    json const& boundary,
    std::unordered_map<std::string, int> const& dof_table,
    double const generate_time_step)
{
    if (std::string const& type = boundary["type"]; type == "traction")
    {
        for (auto it = dof_table.begin(); it != dof_table.end(); ++it)
        {
            if (boundary.find(it->first) != boundary.end())
            {
                auto const dof_offset = it->second;

                for (auto const& mesh : submeshes)
                {
                    boundary_meshes.emplace_back(std::in_place_type_t<traction>{},
                                                 mesh.all_node_indices(),
                                                 2 * mesh.all_node_indices() + dof_offset,
                                                 coordinates,
                                                 boundary,
                                                 it->first,
                                                 generate_time_step,
                                                 mesh.topology(),
                                                 simulation_data);
                }
            }
        }
    }
    else if (type == "body_force")
    {
        for (auto it = dof_table.begin(); it != dof_table.end(); ++it)
        {
            if (boundary.find(it->first) != boundary.end())
            {
                auto const dof_offset = it->second;

                for (auto const& mesh : submeshes)
                {
                    boundary_meshes.emplace_back(std::in_place_type_t<body_force>{},
                                                 mesh.all_node_indices(),
                                                 2 * mesh.all_node_indices() + dof_offset,
                                                 coordinates,
                                                 boundary,
                                                 it->first,
                                                 generate_time_step,
                                                 mesh.topology(),
                                                 simulation_data);
                }
            }
        }
    }
    else
    {
        throw std::domain_error("Need to specify a boundary type \"traction\", "
                                "\"pressure\" or "
                                "\"body_force\"");
    }
}
}
