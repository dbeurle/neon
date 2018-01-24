
#include "NonFollowerLoad.hpp"

#include "interpolations/interpolation_factory.hpp"

#include <utility>

#include "io/json.hpp"

namespace neon::mechanical::solid
{
NonFollowerLoadBoundary::NonFollowerLoadBoundary(
    std::shared_ptr<MaterialCoordinates>& material_coordinates,
    std::vector<Submesh> const& submeshes,
    json const& simulation_data,
    json const& boundary,
    std::unordered_map<std::string, int> const& dof_table)
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
                throw std::runtime_error("x, y or z are acceptable coordinates\n");
            }

            auto const dof_offset = dof_table.find(it.key())->second;

            auto& [is_dof_active, boundary_meshes] = nonfollower_load[dof_offset];

            is_dof_active = true;

            for (auto const& mesh : submeshes)
            {
                boundary_meshes.emplace_back(std::in_place_type_t<traction>{},
                                             make_surface_interpolation(mesh.topology(),
                                                                        simulation_data),
                                             mesh.connectivities(),
                                             filter_dof_list(3, dof_offset, mesh.connectivities()),
                                             material_coordinates,
                                             boundary["Time"],
                                             it.value());
            }
        }
    }
    else if (type == "Pressure")
    {
        auto& [is_dof_active, boundary_meshes] = nonfollower_load[0];

        is_dof_active = true;

        for (auto const& mesh : submeshes)
        {
            boundary_meshes.emplace_back(std::in_place_type_t<pressure>{},
                                         make_surface_interpolation(mesh.topology(), simulation_data),
                                         mesh.connectivities(),
                                         allocate_dof_list(3, mesh.connectivities()),
                                         material_coordinates,
                                         boundary["Time"],
                                         boundary["Values"]);
        }
    }
    else if (type == "BodyForce")
    {
        for (json::const_iterator it = values.begin(); it != values.end(); ++it)
        {
            if (dof_table.find(it.key()) == dof_table.end())
            {
                throw std::runtime_error("x, y or z are acceptable coordinates\n");
            }
            auto const dof_offset = dof_table.find(it.key())->second;

            auto& [is_dof_active, boundary_meshes] = nonfollower_load[dof_offset];

            is_dof_active = true;

            for (auto const& mesh : submeshes)
            {
                boundary_meshes.emplace_back(std::in_place_type_t<body_force>{},
                                             make_volume_interpolation(mesh.topology(),
                                                                       simulation_data),
                                             mesh.connectivities(),
                                             filter_dof_list(3, dof_offset, mesh.connectivities()),
                                             material_coordinates,
                                             boundary["Time"],
                                             it.value());
            }
        }
    }
    else
    {
        throw std::runtime_error("Need to specify a boundary type \"Traction\", \"Pressure\" or "
                                 "\"BodyForce\"");
    }
}
}
