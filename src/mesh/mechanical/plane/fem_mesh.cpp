
#include "mesh/mechanical/plane/fem_mesh.hpp"

#include "mesh/basic_mesh.hpp"
#include "mesh/unique_dof_allocator.hpp"
#include "io/json.hpp"

#include <chrono>
#include <exception>
#include <memory>
#include <numeric>

#include <termcolor/termcolor.hpp>

namespace neon::mechanical::plane
{
fem_mesh::fem_mesh(basic_mesh const& basic_mesh,
                   json const& material_data,
                   json const& simulation_data,
                   double const generate_time_step)
    : coordinates(std::make_shared<material_coordinates>(basic_mesh.coordinates())),
      internal_forces{coordinates->size() * traits::dofs_per_node},
      generate_time_step{generate_time_step}
{
    check_boundary_conditions(simulation_data["BoundaryConditions"]);

    for (auto const& submesh : basic_mesh.meshes(simulation_data["Name"]))
    {
        submeshes.emplace_back(material_data, simulation_data, coordinates, submesh);
    }
    allocate_boundary_conditions(simulation_data, basic_mesh);
}

bool fem_mesh::is_symmetric() const
{
    return std::all_of(begin(submeshes), end(submeshes), [](auto const& submesh) {
        return submesh.constitutive().is_symmetric();
    });
}

void fem_mesh::update_internal_variables(vector const& u, double const time_step_size)
{
    auto const start = std::chrono::high_resolution_clock::now();

    coordinates->update_current_xy_configuration(u);

    for (auto& submesh : submeshes) submesh.update_internal_variables(time_step_size);

    auto const end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;

    std::cout << std::string(6, ' ') << "Internal variable update took " << elapsed_seconds.count()
              << "s\n";
}

void fem_mesh::save_internal_variables(bool const have_converged)
{
    for (auto& submesh : submeshes) submesh.save_internal_variables(have_converged);
}

bool fem_mesh::is_nonfollower_load(std::string const& boundary_type) const
{
    return boundary_type == "Traction" || boundary_type == "Pressure" || boundary_type == "BodyForce";
}

void fem_mesh::allocate_boundary_conditions(json const& simulation_data, basic_mesh const& basic_mesh)
{
    // Populate the boundary conditions and their corresponding mesh
    for (auto const& boundary : simulation_data["BoundaryConditions"])
    {
        auto const& boundary_name = boundary["Name"].get<std::string>();
        auto const& boundary_type = boundary["Type"].get<std::string>();

        if (boundary_type == "Displacement")
        {
            this->allocate_displacement_boundary(boundary, basic_mesh);
        }
        else if (is_nonfollower_load(boundary_type))
        {
            nonfollower_loads.emplace(boundary_name,
                                      nonfollower_load_boundary(coordinates,
                                                                basic_mesh.meshes(boundary_name),
                                                                simulation_data,
                                                                boundary,
                                                                dof_table,
                                                                generate_time_step));
        }
        else
        {
            throw std::domain_error("BoundaryCondition \"" + boundary_type + "\" is not recognised");
        }
    }
}

void fem_mesh::allocate_displacement_boundary(json const& boundary, basic_mesh const& basic_mesh)
{
    std::string const& boundary_name = boundary["Name"];

    for (auto const& [dof_key, dof_offset] : dof_table)
    {
        if (boundary.count(dof_key))
        {
            auto boundary_dofs = unique_dof_allocator<traits::dofs_per_node>(
                basic_mesh.meshes(boundary_name));

            // Offset the degrees of freedom on the boundary
            std::transform(begin(boundary_dofs),
                           end(boundary_dofs),
                           begin(boundary_dofs),
                           [&](auto const dof) { return dof + dof_offset; });

            displacement_bcs[boundary_name].emplace_back(boundary_dofs,
                                                         boundary,
                                                         dof_key,
                                                         generate_time_step);
        }
    }
}

std::vector<double> fem_mesh::time_history() const
{
    std::set<double> history;

    // Append time history from each boundary condition
    for (auto const& [key, boundaries] : displacement_bcs)
    {
        for (auto const& boundary : boundaries)
        {
            auto const& times = boundary.time_history();

            history.insert(begin(times), end(times));
        }
    }
    for (auto const& [key, nonfollower_load] : nonfollower_loads)
    {
        for (auto const& [is_dof_active, boundaries] : nonfollower_load.interface())
        {
            if (!is_dof_active) continue;

            for (auto const& boundary_variant : boundaries)
            {
                std::visit(
                    [&](auto const& surface_mesh) {
                        auto const& times = surface_mesh.time_history();

                        history.insert(begin(times), end(times));
                    },
                    boundary_variant);
            }
        }
    }
    return {begin(history), end(history)};
}

void fem_mesh::check_boundary_conditions(json const& boundary_data) const
{
    for (auto const& boundary : boundary_data)
    {
        for (auto const& mandatory_field : {"Name", "Time", "Type"})
        {
            if (!boundary.count(mandatory_field))
            {
                throw std::domain_error("\"" + std::string(mandatory_field)
                                        + "\" was not specified in \"BoundaryCondition\".");
            }
        }
    }
}
}
