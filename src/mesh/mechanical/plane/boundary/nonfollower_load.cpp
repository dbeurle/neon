
#include "nonfollower_load.hpp"

#include "geometry/Projection.hpp"
#include "math/transform_expand.hpp"
#include "io/json.hpp"

#include <Eigen/Geometry>

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

    auto const& values = boundary["Values"];

    if (auto const& type = boundary["Type"].get<std::string>(); type == "Traction")
    {
        for (json::const_iterator it = values.begin(); it != values.end(); ++it)
        {
            if (dof_table.find(it.key()) == dof_table.end())
            {
                throw std::domain_error("x, y or z are acceptable coordinates\n");
            }

            auto const dof_offset = dof_table.find(it.key())->second;

            auto& [is_dof_active, boundary_meshes] = nonfollower_load[dof_offset];

            is_dof_active = true;

            for (auto const& mesh : submeshes)
            {
                boundary_meshes.emplace_back(std::in_place_type_t<traction>{},
                                             make_line_interpolation(mesh.topology(), simulation_data),
                                             mesh.element_connectivity(),
                                             2 * mesh.element_connectivity() + dof_offset,
                                             material_coordinates,
                                             boundary,
                                             it.key(),
                                             generate_time_step);
            }
        }
    }
    else if (type == "BodyForce")
    {
        for (json::const_iterator it = values.begin(); it != values.end(); ++it)
        {
            if (dof_table.find(it.key()) == dof_table.end())
            {
                throw std::domain_error("x or y are acceptable coordinates\n");
            }
            auto const dof_offset = dof_table.find(it.key())->second;

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
                                             it.key(),
                                             generate_time_step);
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
