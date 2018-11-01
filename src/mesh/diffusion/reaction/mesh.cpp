
#include "mesh.hpp"

#include "mesh/basic_mesh.hpp"
#include "mesh/dof_allocator.hpp"
#include "io/json.hpp"
#include "io/post/variable_string_adapter.hpp"
#include "io/post/node_averaged_variables.hpp"

#include <termcolor/termcolor.hpp>

#include <chrono>
#include <exception>
#include <numeric>

static bool is_nodal_variable(std::string const& name) { return name == "specie"; }

namespace neon::diffusion::reaction
{
mesh::mesh(basic_mesh const& basic_mesh, json const& material_data, json const& mesh_data)
    : coordinates(std::make_shared<material_coordinates>(basic_mesh.coordinates())),
      writer(std::make_unique<io::vtk_file_output>(mesh_data["name"], mesh_data["visualisation"]))
{
    check_boundary_conditions(mesh_data["boundaries"]);

    std::string const& simulation_name = mesh_data["name"];

    writer->coordinates(coordinates->coordinates());

    for (auto const& submesh : basic_mesh.meshes(simulation_name))
    {
        submeshes.emplace_back(material_data, mesh_data, coordinates, submesh);

        writer->mesh(submesh.all_node_indices(), submesh.topology());
    }
    allocate_boundary_conditions(mesh_data, basic_mesh);

    allocate_variable_names();
}

void mesh::update_internal_variables(vector const& u, double const time_step_size)
{
    auto const start = std::chrono::steady_clock::now();

    temperature = u;

    for (auto& submesh : submeshes)
    {
        submesh.update_internal_variables(time_step_size);
    }

    auto const end = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;

    std::cout << std::string(6, ' ') << "Internal variable update took " << elapsed_seconds.count()
              << "s\n";
}

void mesh::save_internal_variables(bool const have_converged)
{
    for (auto& submesh : submeshes)
    {
        submesh.save_internal_variables(have_converged);
    }
}

void mesh::allocate_boundary_conditions(json const& mesh_data, basic_mesh const& basic_mesh)
{
    auto const& boundary_data = mesh_data["boundaries"];

    // Populate the boundary meshes
    for (auto const& boundary : boundary_data)
    {
        std::string const& boundary_name = boundary["name"];

        if (boundary.find("time") == boundary.end())
        {
            throw std::domain_error("boundaries requires a \"time\" field.");
        }

        if (std::string const& boundary_type = boundary["type"]; boundary_type == "concentration")
        {
            if (boundary.find("value") == boundary.end())
            {
                throw std::domain_error("boundaries \"" + boundary_type
                                        + "\" requires a \"value\" field.");
            }

            dirichlet_bcs[boundary_name].emplace_back(unique_dof_allocator<1>(
                                                          basic_mesh.meshes(boundary_name)),
                                                      boundary["time"],
                                                      boundary["value"]);
        }
        else if (boundary_type == "flux")
        {
            if (boundary_type == "flux" && boundary.find("value") == boundary.end())
            {
                throw std::domain_error("boundaries \"" + boundary_type
                                        + "\" requires a \"value\" field.");
            }
            boundary_meshes[boundary_name].emplace_back(coordinates,
                                                        basic_mesh.meshes(boundary_name),
                                                        boundary,
                                                        mesh_data);
        }
        else
        {
            throw std::domain_error("boundaries \"" + boundary_type
                                    + "\" is not recognised.  For reaction-diffusion simulations "
                                      "\"concentration\" and \"flux\" are valid.");
        }
    }
}

void mesh::check_boundary_conditions(json const& boundary_data) const
{
    for (auto const& boundary : boundary_data)
    {
        for (auto const& required_field : {"name", "type"})
        {
            if (boundary.find(required_field) == boundary.end())
            {
                throw std::domain_error("Missing " + std::string(required_field) + " in boundaries\n");
            }
        }
    }
}

void mesh::write(std::int32_t const time_step, double const current_time)
{
    // nodal variables
    if (writer->is_output_requested("temperature"))
    {
        writer->field("temperature", temperature, 1);
    }
    // internal variables extrapolated to the nodes
    for (auto const& output_variable : output_variables)
    {
        std::visit(
            [this](auto&& output) {
                using T = std::decay_t<decltype(output)>;
                if constexpr (std::is_same_v<T, variable::scalar>)
                {
                    writer->field(convert(output),
                                  average_internal_variable(submeshes,
                                                            coordinates->size(),
                                                            convert(output),
                                                            output),
                                  1);
                }
                else if constexpr (std::is_same_v<T, variable::second>)
                {
                    writer->field(convert(output),
                                  average_internal_variable(submeshes,
                                                            coordinates->size() * 9,
                                                            convert(output),
                                                            output),
                                  9);
                }
            },
            output_variable);
    }
    writer->write(time_step, current_time);
}

void mesh::allocate_variable_names()
{
    for (auto const& name : writer->outputs())
    {
        if (is_nodal_variable(name))
        {
            continue;
        }
        output_variables.emplace_back(variable::convert(name));
    }

    // check if output variables exist in all the submeshes
    for (auto const& output_variable : output_variables)
    {
        std::visit(
            [this](auto&& output) {
                using T = std::decay_t<decltype(output)>;
                if constexpr (std::is_same_v<T, variable::scalar>)
                {
                    if (std::none_of(begin(submeshes), end(submeshes), [&output](auto const& submesh) {
                            return submesh.internal_variables().has(output);
                        }))
                    {
                        throw std::domain_error("Internal variables do not have the requested "
                                                "scalar variable ("
                                                + std::to_string(static_cast<short>(output)) + ")");
                    }
                }
                else if constexpr (std::is_same_v<T, variable::second>)
                {
                    if (std::none_of(begin(submeshes), end(submeshes), [&output](auto const& submesh) {
                            return submesh.internal_variables().has(output);
                        }))
                    {
                        throw std::domain_error("Internal variables do not have the requested "
                                                "second order tensor variable ("
                                                + std::to_string(static_cast<short>(output)) + ")");
                    }
                }
            },
            output_variable);
    }
}
}
