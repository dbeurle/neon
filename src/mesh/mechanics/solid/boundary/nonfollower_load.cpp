
#include "nonfollower_load.hpp"

#include "interpolations/interpolation_factory.hpp"
#include "io/json.hpp"
#include "traits/mechanics.hpp"
#include "mesh/dof_allocator.hpp"

#include <utility>

namespace neon::mechanics::solid
{
nonfollower_load_boundary::nonfollower_load_boundary(
    std::shared_ptr<material_coordinates>& material_coordinates,
    std::vector<basic_submesh> const& submeshes,
    json const& simulation_data,
    json const& boundary_data,
    std::unordered_map<std::string, int> const& dof_table,
    double const generate_time_step)
{
    if (std::string const& type = boundary_data["type"]; type == "traction")
    {
        for (auto const& [dof_name, dof_offset] : dof_table)
        {
            if (boundary_data.find(dof_name) != end(boundary_data))
            {
                for (auto const& submesh : submeshes)
                {
                    boundary_meshes.emplace_back(std::in_place_type_t<traction>{},
                                                 submesh.all_node_indices(),
                                                 3 * submesh.all_node_indices() + dof_offset,
                                                 material_coordinates,
                                                 boundary_data,
                                                 dof_name,
                                                 generate_time_step,
                                                 submesh.topology(),
                                                 simulation_data);
                }
            }
        }
    }
    else if (type == "pressure")
    {
        if (boundary_data.find("value") == end(boundary_data))
        {
            throw std::domain_error("Pressure boundary condition must specify a \"value\" array");
        }

        for (auto const& submesh : submeshes)
        {
            auto const& node_indices = submesh.all_node_indices();

            indices dof_indices(3 * node_indices.rows(), node_indices.cols());

            dof_allocator(node_indices,
                          dof_indices,
                          traits<theory::solid, discretisation::linear, true>::dofs_per_node);

            boundary_meshes.emplace_back(std::in_place_type_t<pressure>{},
                                         node_indices,
                                         dof_indices,
                                         material_coordinates,
                                         boundary_data["time"],
                                         boundary_data["value"],
                                         submesh.topology(),
                                         simulation_data);
        }
    }
    else if (type == "body_force")
    {
        for (auto const& [dof_name, dof_offset] : dof_table)
        {
            if (boundary_data.find(dof_name) != end(boundary_data))
            {
                for (auto const& submesh : submeshes)
                {
                    boundary_meshes.emplace_back(std::in_place_type_t<body_force>{},
                                                 submesh.all_node_indices(),
                                                 3 * submesh.all_node_indices() + dof_offset,
                                                 material_coordinates,
                                                 boundary_data,
                                                 dof_name,
                                                 generate_time_step,
                                                 submesh.topology(),
                                                 simulation_data);
                }
            }
        }
    }
    else if (type == "nodal_force")
    {
        for (auto const& [dof_name, dof_offset] : dof_table)
        {
            if (boundary_data.find(dof_name) != end(boundary_data))
            {
                for (auto const& submesh : submeshes)
                {
                    // create linear dof indices
                    auto node_indices = submesh.unique_node_indices();

                    // Offset the degrees of freedom on the boundary
                    std::transform(begin(node_indices),
                                   end(node_indices),
                                   begin(node_indices),
                                   [&, dof_offset = std::ref(dof_offset)](auto const node) {
                                       return node * 3 + dof_offset;
                                   });

                    nodal_values.emplace_back(node_indices, boundary_data, dof_name, generate_time_step);
                }
            }
        }
    }
    else
    {
        throw std::domain_error("Need to specify a boundary type \"traction\", \"pressure\" or "
                                "\"body_force\"");
    }
}
}
