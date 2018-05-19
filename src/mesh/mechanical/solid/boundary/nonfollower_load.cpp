
#include "nonfollower_load.hpp"

#include "interpolations/interpolation_factory.hpp"
#include "io/json.hpp"
#include "traits/mechanics.hpp"
#include "math/transform_expand.hpp"

#include <utility>

namespace neon::mechanical::solid
{
nonfollower_load_boundary::nonfollower_load_boundary(
    std::shared_ptr<material_coordinates>& material_coordinates,
    std::vector<basic_submesh> const& submeshes,
    json const& simulation_data,
    json const& boundary_data,
    std::unordered_map<std::string, int> const& dof_table,
    double const generate_time_step)
{
    if (std::string const& type = boundary_data["Type"]; type == "Traction")
    {
        for (auto const& [dof_name, dof_offset] : dof_table)
        {
            if (boundary_data.count(dof_name))
            {
                for (auto const& mesh : submeshes)
                {
                    boundary_meshes.emplace_back(std::in_place_type_t<traction>{},
                                                 make_surface_interpolation(mesh.topology(),
                                                                            simulation_data),
                                                 mesh.all_node_indices(),
                                                 3 * mesh.all_node_indices() + dof_offset,
                                                 material_coordinates,
                                                 boundary_data,
                                                 dof_name,
                                                 generate_time_step);
                }
            }
        }
    }
    else if (type == "Pressure")
    {
        for (auto const& mesh : submeshes)
        {
            auto const& node_indices = mesh.all_node_indices();

            indices dof_indices(3 * node_indices.rows(), node_indices.cols());

            for (indices::Index i{0}; i < node_indices.cols(); ++i)
            {
                transform_expand_view(node_indices(Eigen::placeholders::all, i),
                                      dof_indices(Eigen::placeholders::all, i),
                                      traits<theory::solid, discretisation::linear, true>::dof_order);
            }

            boundary_meshes.emplace_back(std::in_place_type_t<pressure>{},
                                         make_surface_interpolation(mesh.topology(), simulation_data),
                                         node_indices,
                                         dof_indices,
                                         material_coordinates,
                                         boundary_data["Time"],
                                         boundary_data["Value"]);
        }
    }
    else if (type == "BodyForce")
    {
        for (auto const& [dof_name, dof_offset] : dof_table)
        {
            if (boundary_data.count(dof_name))
            {
                for (auto const& mesh : submeshes)
                {
                    boundary_meshes.emplace_back(std::in_place_type_t<body_force>{},
                                                 make_volume_interpolation(mesh.topology(),
                                                                           simulation_data),
                                                 mesh.all_node_indices(),
                                                 3 * mesh.all_node_indices() + dof_offset,
                                                 material_coordinates,
                                                 boundary_data,
                                                 dof_name,
                                                 generate_time_step);
                }
            }
        }
    }
    else if (type == "NodalForce")
    {
        for (auto const& [dof_name, dof_offset] : dof_table)
        {
            if (boundary_data.count(dof_name))
            {
                for (auto const& mesh : submeshes)
                {
                    // create linear dof indices
                    auto node_indices = mesh.unique_node_indices();

                    // Offset the degrees of freedom on the boundary
                    std::transform(begin(node_indices),
                                   end(node_indices),
                                   begin(node_indices),
                                   [&](auto const node) { return node * 3 + dof_offset; });

                    nodal_values.emplace_back(node_indices, boundary_data, dof_name, generate_time_step);
                }
            }
        }
    }
    else
    {
        throw std::domain_error("Need to specify a boundary type \"Traction\", \"Pressure\" or "
                                "\"BodyForce\"");
    }
}
}
