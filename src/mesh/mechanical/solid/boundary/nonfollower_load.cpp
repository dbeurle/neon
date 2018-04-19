
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
                                                 make_surface_interpolation(mesh.topology(),
                                                                            simulation_data),
                                                 mesh.element_connectivity(),
                                                 3 * mesh.element_connectivity() + dof_offset,
                                                 material_coordinates,
                                                 boundary,
                                                 it->first,
                                                 generate_time_step);
                }
            }
        }
    }
    else if (type == "Pressure")
    {
        auto& [is_dof_active, boundary_meshes] = nonfollower_load[0];

        is_dof_active = true;

        for (auto const& mesh : submeshes)
        {
            auto const& connectivity = mesh.element_connectivity();

            indices dof_list(3 * connectivity.rows(), connectivity.cols());

            for (indices::Index i{0}; i < connectivity.cols(); ++i)
            {
                transform_expand_view(connectivity(Eigen::placeholders::all, i),
                                      dof_list(Eigen::placeholders::all, i),
                                      traits<theory::solid, discretisation::linear, true>::dof_order);
            }

            boundary_meshes.emplace_back(std::in_place_type_t<pressure>{},
                                         make_surface_interpolation(mesh.topology(), simulation_data),
                                         mesh.element_connectivity(),
                                         dof_list,
                                         material_coordinates,
                                         boundary["Time"],
                                         boundary["Value"]);
        }
    }
    else if (type == "BodyForce")
    {
        for (auto it = dof_table.begin(); it != dof_table.end(); ++it)
        {
            if (boundary.count(it->first))
            {
                auto const dof_offset = it->second;

                auto& [is_dof_active, boundary_meshes] = nonfollower_load[dof_offset];

                is_dof_active = true;

                for (auto const& mesh : submeshes)
                {
                    boundary_meshes.emplace_back(std::in_place_type_t<body_force>{},
                                                 make_volume_interpolation(mesh.topology(),
                                                                           simulation_data),
                                                 mesh.element_connectivity(),
                                                 3 * mesh.element_connectivity() + dof_offset,
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
        throw std::domain_error("Need to specify a boundary type \"Traction\", \"Pressure\" or "
                                "\"BodyForce\"");
    }
}
}
