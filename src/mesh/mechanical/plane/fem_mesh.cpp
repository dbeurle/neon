
#include "mesh/mechanical/plane/fem_mesh.hpp"

#include "mesh/basic_mesh.hpp"
#include "mesh/mesh_dof_filter.hpp"

#include <chrono>
#include <exception>
#include <memory>
#include <numeric>

#include "io/json.hpp"
#include <termcolor/termcolor.hpp>

#include <range/v3/action/join.hpp>
#include <range/v3/action/sort.hpp>
#include <range/v3/action/transform.hpp>
#include <range/v3/action/unique.hpp>
#include <range/v3/view/transform.hpp>

namespace neon::mechanical::plane
{
fem_mesh::fem_mesh(basic_mesh const& basic_mesh,
                   json const& material_data,
                   json const& simulation_data,
                   double const generate_time_step)
    : mesh_coordinates(std::make_shared<material_coordinates>(basic_mesh.coordinates())),
      generate_time_step{generate_time_step}
{
    check_boundary_conditions(simulation_data["BoundaryConditions"]);

    auto const& simulation_name = simulation_data["Name"].get<std::string>();

    for (auto const& submesh : basic_mesh.meshes(simulation_name))
    {
        submeshes.emplace_back(material_data, simulation_data, mesh_coordinates, submesh);
    }
    allocate_boundary_conditions(simulation_data, basic_mesh);
}

bool fem_mesh::is_symmetric() const
{
    for (auto const& submesh : submeshes)
    {
        if (!submesh.constitutive().is_symmetric()) return false;
    }
    return true;
}

void fem_mesh::update_internal_variables(vector const& u, double const time_step_size)
{
    auto const start = std::chrono::high_resolution_clock::now();

    mesh_coordinates->update_current_xy_configuration(u);

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
    auto const& boundary_data = simulation_data["BoundaryConditions"];

    // Populate the boundary conditions and their corresponding mesh
    for (auto const& boundary : boundary_data)
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
                                      nonfollower_load_boundary(mesh_coordinates,
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
    using namespace ranges;

    auto const& boundary_name = boundary["Name"].get<std::string>();

    auto const dirichlet_dofs = mesh_dof_filter<2>(basic_mesh.meshes(boundary_name));

    for (auto it = dof_table.begin(); it != dof_table.end(); ++it)
    {
        if (boundary.count(it->first))
        {
            auto const dof_offset = it->second;

            // Offset the degrees of freedom on the boundary
            auto const boundary_dofs = view::transform(dirichlet_dofs, [&](auto const& dof) {
                return dof + dof_offset;
            });

            displacement_bcs[boundary_name].emplace_back(boundary_dofs,
                                                         boundary,
                                                         it->first,
                                                         generate_time_step);
        }
    }
}

std::vector<double> fem_mesh::time_history() const
{
    std::vector<double> history;

    // Append time history from each boundary condition
    for (auto const& [key, boundaries] : displacement_bcs)
    {
        for (auto const& boundary : boundaries)
        {
            history |= ranges::action::push_back(boundary.time_history());
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
                        history |= ranges::action::push_back(surface_mesh.time_history());
                    },
                    boundary_variant);
            }
        }
    }
    return std::move(history) | ranges::action::sort | ranges::action::unique;
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
