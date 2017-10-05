
#include "NonFollowerLoad.hpp"

#include <iostream>
#include <utility>

#include <json/value.h>

namespace neon::solid
{
Pressure::Pressure(std::unique_ptr<SurfaceInterpolation>&& sf,
                   std::vector<List> const& nodal_connectivity,
                   std::shared_ptr<MaterialCoordinates>& material_coordinates,
                   Json::Value const& time_history,
                   Json::Value const& load_history,
                   int const nodal_dofs)
    : Traction(std::move(sf),
               nodal_connectivity,
               material_coordinates,
               time_history,
               load_history,
               0,
               nodal_dofs)
{
    // The inclusion of pressure results in at most three tractions
    dof_list = allocate_dof_list(3, nodal_connectivity);
}

std::tuple<List const&, Vector> Pressure::external_force(int const element,
                                                         double const load_factor) const
{
    auto const X = material_coordinates->initial_configuration(nodal_connectivity[element]);

    auto const X_surface = sf->project_to_plane(X);

    auto const p = interpolate_prescribed_load(load_factor);

    // Perform the computation of the external load vector
    RowMatrix f_ext = -p
                      * sf->quadrature()
                            .integrate(RowMatrix::Zero(X.cols(), 3).eval(),
                                       [&](auto const& femval, auto const& l) -> RowMatrix {
                                           auto const & [ N, dN ] = femval;

                                           auto const j = (X_surface * dN).determinant();

                                           Vector3 const x_xi = (X * dN).col(0);
                                           Vector3 const x_eta = (X * dN).col(1);

                                           Vector3 const normal = x_xi.cross(x_eta).normalized();

                                           return N * normal.transpose() * j;
                                       });

    // Map the matrix back to a vector for the assembly operator
    return {dof_list[element], Eigen::Map<RowMatrix>(f_ext.data(), X.cols() * 3, 1)};
}

NonFollowerLoadBoundary::NonFollowerLoadBoundary(
    std::shared_ptr<MaterialCoordinates>& material_coordinates,
    std::vector<SubMesh> const& submeshes,
    Json::Value const& simulation_data,
    Json::Value const& boundary,
    std::unordered_map<std::string, int> const& dof_table)
{
    auto const& type = boundary["Type"].asString();

    for (auto & [ is_dof_active, var ] : nonfollower_load)
    {
        is_dof_active = false;
    }

    if (type == "Traction")
    {
        for (auto const& name : boundary["Values"].getMemberNames())
        {
            if (dof_table.find(name) == dof_table.end())
            {
                throw std::runtime_error("x, y or z are acceptable coordinates\n");
            }

            auto const dof_offset = dof_table.find(name)->second;

            auto & [ is_dof_active, boundary_meshes ] = nonfollower_load[dof_offset];

            is_dof_active = true;

            for (auto const& mesh : submeshes)
            {
                boundary_meshes.emplace_back(std::in_place_type_t<Traction>{},
                                             make_surface_interpolation(mesh.topology(),
                                                                        simulation_data),
                                             mesh.connectivities(),
                                             material_coordinates,
                                             boundary["Time"],
                                             boundary["Values"][name],
                                             dof_offset,
                                             3);
            }
        }
    }
    else if (type == "Pressure")
    {
        auto & [ is_dof_active, boundary_meshes ] = nonfollower_load[0];

        is_dof_active = true;

        for (auto const& mesh : submeshes)
        {
            boundary_meshes.emplace_back(std::in_place_type_t<Pressure>{},
                                         make_surface_interpolation(mesh.topology(), simulation_data),
                                         mesh.connectivities(),
                                         material_coordinates,
                                         boundary["Time"],
                                         boundary["Values"],
                                         3);
        }
    }
    else if (type == "BodyForce")
    {
        for (auto const& name : boundary["Values"].getMemberNames())
        {
            if (dof_table.find(name) == dof_table.end())
            {
                throw std::runtime_error("x, y or z are acceptable coordinates\n");
            }
            auto const dof_offset = dof_table.find(name)->second;

            auto & [ is_dof_active, boundary_meshes ] = nonfollower_load[dof_offset];

            is_dof_active = true;

            for (auto const& mesh : submeshes)
            {
                boundary_meshes.emplace_back(std::in_place_type_t<BodyForce>{},
                                             make_volume_interpolation(mesh.topology(),
                                                                       simulation_data),
                                             mesh.connectivities(),
                                             material_coordinates,
                                             boundary["Time"],
                                             boundary["Values"][name],
                                             dof_offset,
                                             3);
            }
        }
    }
}
}
