
#include "nonfollower_load.hpp"

#include "geometry/Projection.hpp"
#include "math/transform_expand.hpp"
#include "io/json.hpp"

namespace neon::mechanical::plane
{
nonfollower_load_boundary::nonfollower_load_boundary(
    std::shared_ptr<material_coordinates>& material_coordinates,
    std::vector<basic_submesh> const& submeshes,
    json const& simulation_data,
    json const& boundary,
    std::unordered_map<std::string, int> const& dof_table,
    double const generate_time_step)
{
    for (auto& [is_dof_active, var] : nonfollower_load)
    {
        is_dof_active = false;
    }

    if (auto const& type = boundary["Type"].get<std::string>(); type == "Traction")
    {
        for (auto it = dof_table.begin(); it != dof_table.end(); ++it)
        {
            if (boundary.count(it->first))
            {
                auto const& dof_offset = it->second;

                auto& [is_dof_active, boundary_meshes] = nonfollower_load[dof_offset];

                is_dof_active = true;

                for (auto const& mesh : submeshes)
                {
                    boundary_meshes.emplace_back(std::in_place_type_t<traction>{},
                                                 make_line_interpolation(mesh.topology(),
                                                                         simulation_data),
                                                 mesh.element_connectivity(),
                                                 2 * mesh.element_connectivity() + dof_offset,
                                                 material_coordinates,
                                                 boundary,
                                                 it->first,
                                                 generate_time_step);
                }
            }
        }
    }
    else if (type == "BodyForce")
    {
        for (auto it = dof_table.begin(); it != dof_table.end(); ++it)
        {
            if (boundary.count(it->first))
            {
                auto const& dof_offset = it->second;

                auto& [is_dof_active, boundary_meshes] = nonfollower_load[dof_offset];

                is_dof_active = true;

                for (auto const& mesh : submeshes)
                {
                    boundary_meshes.emplace_back(std::in_place_type_t<body_force>{},
                                                 make_surface_interpolation(mesh.topology(),
                                                                            simulation_data),
                                                 mesh.element_connectivity(),
                                                 2 * mesh.element_connectivity() + dof_offset,
                                                 material_coordinates,
                                                 boundary,
                                                 it->first,
                                                 generate_time_step);
                }
            }
        }
    }
    else
    {
        throw std::domain_error("Need to specify a boundary type \"Traction\", "
                                "\"Pressure\" or "
                                "\"BodyForce\"");
    }
}
}
