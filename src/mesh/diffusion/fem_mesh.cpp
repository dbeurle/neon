
#include "fem_mesh.hpp"

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
#include <range/v3/action/unique.hpp>
#include <range/v3/view/transform.hpp>

namespace neon::diffusion
{
fem_mesh::fem_mesh(basic_mesh const& basic_mesh, json const& material_data, json const& mesh_data)
    : mesh_coordinates(std::make_shared<material_coordinates>(basic_mesh.coordinates()))
{
    check_boundary_conditions(mesh_data["BoundaryConditions"]);

    auto const& simulation_name = mesh_data["Name"].get<std::string>();

    for (auto const& submesh : basic_mesh.meshes(simulation_name))
    {
        submeshes.emplace_back(material_data, mesh_data, mesh_coordinates, submesh);
    }
    allocate_boundary_conditions(mesh_data, basic_mesh);
}

void fem_mesh::update_internal_variables(vector const& u, double const time_step_size)
{
    auto const start = std::chrono::high_resolution_clock::now();

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

void fem_mesh::allocate_boundary_conditions(json const& mesh_data, basic_mesh const& basic_mesh)
{
    auto const& boundary_data = mesh_data["BoundaryConditions"];

    // Populate the boundary meshes
    for (auto const& boundary : boundary_data)
    {
        std::string const& boundary_name = boundary["Name"];

        if (!boundary.count("Time"))
        {
            throw std::domain_error("BoundaryCondition requires a \"Time\" field.");
        }

        if (auto const& boundary_type = boundary["Type"].get<std::string>();
            boundary_type == "Temperature")
        {
            if (!boundary.count("Value"))
            {
                throw std::domain_error("BoundaryCondition \"" + boundary_type
                                        + "\" requires a \"Value\" field.");
            }

            dirichlet_bcs[boundary_name].emplace_back(mesh_dof_filter<1>(
                                                          basic_mesh.meshes(boundary_name)),
                                                      boundary["Time"],
                                                      boundary["Value"]);
        }
        else if (boundary_type == "HeatFlux" || boundary_type == "NewtonCooling")
        {
            if (boundary_type == "HeatFlux" && !boundary.count("Value"))
            {
                throw std::domain_error("BoundaryCondition \"" + boundary_type
                                        + "\" requires a \"Value\" field.");
            }
            else if (boundary_type == "NewtonCooling"
                     && (!boundary.count("HeatTransferCoefficient")
                         || !boundary.count("AmbientTemperature")))
            {
                throw std::domain_error("BoundaryCondition \"" + boundary_type
                                        + "\" requires a \"HeatTransferCoefficient\" and "
                                          "\"AmbientTemperature\" "
                                          "field.");
            }
            boundary_meshes[boundary_name].emplace_back(mesh_coordinates,
                                                        basic_mesh.meshes(boundary_name),
                                                        boundary,
                                                        mesh_data);
        }
        else
        {
            throw std::domain_error("BoundaryCondition \"" + boundary_type
                                    + "\" is not recognised.  For thermal simulations "
                                      "\"Temperature\", \"HeatFlux\" "
                                      "and \"NewtonCooling\" are valid.");
        }
    }
}

void fem_mesh::check_boundary_conditions(json const& boundary_data) const
{
    for (auto const& boundary : boundary_data)
    {
        for (auto const& required_field : {"Name", "Type"})
        {
            if (!boundary.count(required_field))
            {
                throw std::domain_error("Missing " + std::string(required_field)
                                        + " in BoundaryConditions\n");
            }
        }
    }
}
}
