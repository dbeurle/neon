
#include "fem_mesh.hpp"

#include "mesh/basic_mesh.hpp"

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

namespace neon::mechanical::solid
{
fem_mesh::fem_mesh(basic_mesh const& basic_mesh, json const& material_data, json const& simulation_data)
    : mesh_coordinates(std::make_shared<material_coordinates>(basic_mesh.coordinates()))
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
    auto start = std::chrono::high_resolution_clock::now();

    mesh_coordinates->update_current_configuration(u);

    for (auto& submesh : submeshes) submesh.update_internal_variables(time_step_size);

    auto end = std::chrono::high_resolution_clock::now();
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
                                      NonFollowerLoadBoundary(mesh_coordinates,
                                                              basic_mesh.meshes(boundary_name),
                                                              simulation_data,
                                                              boundary,
                                                              dof_table));
        }
        else
        {
            throw std::runtime_error("BoundaryCondition \"" + boundary_type + "\" is not recognised");
        }
    }
}

void fem_mesh::allocate_displacement_boundary(json const& boundary, basic_mesh const& basic_mesh)
{
    using namespace ranges;

    auto const& boundary_name = boundary["Name"].get<std::string>();

    auto const dirichlet_dofs = this->filter_dof_list(basic_mesh.meshes(boundary_name));

    auto const& values = boundary["Values"];

    for (json::const_iterator it = values.begin(); it != values.end(); ++it)
    {
        auto const& dof_offset = dof_table.find(it.key())->second;

        if (dof_table.find(it.key()) == dof_table.end())
        {
            throw std::runtime_error("x, y or z are acceptable coordinates\n");
        }

        // Offset the degrees of freedom on the boundary
        auto const boundary_dofs = view::transform(dirichlet_dofs, [&](auto const& dof) {
            return dof + dof_offset;
        });

        displacement_bcs[boundary_name].emplace_back(boundary_dofs, boundary["Time"], it.value());
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
        for (auto const& mandatory_field : {"Name", "Time", "Type", "Values"})
        {
            if (!boundary.count(mandatory_field))
            {
                throw std::runtime_error("\"" + std::string(mandatory_field)
                                         + "\" was not specified in \"BoundaryCondition\".");
            }
        }
    }
}

local_indices fem_mesh::filter_dof_list(std::vector<basic_submesh> const& boundary_mesh) const
{
    using namespace ranges;

    return view::transform(boundary_mesh,
                           [](auto const& submesh) { return submesh.connectivities(); })
           | action::join | action::join | action::sort | action::unique
           | action::transform([=](auto const& i) { return i * 3; });
}
}
