
#include "mesh/mechanical/plane/fem_mesh.hpp"

#include "mesh/basic_mesh.hpp"
#include "mesh/dof_allocator.hpp"
#include "io/json.hpp"
#include "io/post/variable_string_adapter.hpp"
#include "io/post/node_averaged_variables.hpp"

#include <chrono>
#include <exception>
#include <numeric>

#include <termcolor/termcolor.hpp>

static bool is_nodal_variable(std::string const& name)
{
    return name == "displacement" || name == "reaction_force";
}

namespace neon::mechanical::plane
{
fem_mesh::fem_mesh(basic_mesh const& basic_mesh,
                   json const& material_data,
                   json const& simulation_data,
                   double const generate_time_step)
    : coordinates(std::make_shared<material_coordinates>(basic_mesh.coordinates())),
      reaction_forces{coordinates->size() * traits::dofs_per_node},
      generate_time_step{generate_time_step},
      writer(std::make_unique<io::vtk_file_output>(simulation_data["Name"],
                                                   simulation_data["Visualisation"]))
{
    check_boundary_conditions(simulation_data["BoundaryConditions"]);

    writer->coordinates(coordinates->coordinates());

    for (auto const& submesh : basic_mesh.meshes(simulation_data["Name"]))
    {
        submeshes.emplace_back(material_data, simulation_data, coordinates, submesh);

        writer->mesh(submesh.all_node_indices(), submesh.topology());
    }
    allocate_boundary_conditions(simulation_data, basic_mesh);

    allocate_variable_names();
}

bool fem_mesh::is_symmetric() const
{
    return std::all_of(begin(submeshes), end(submeshes), [](auto const& submesh) {
        return submesh.constitutive().is_symmetric();
    });
}

void fem_mesh::update_internal_variables(vector const& u, double const time_step_size)
{
    auto const start = std::chrono::steady_clock::now();

    coordinates->update_current_xy_configuration(u);

    for (auto& submesh : submeshes)
    {
        submesh.update_internal_variables(time_step_size);
    }

    auto const end = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;

    std::cout << std::string(6, ' ') << "Internal variable update took " << elapsed_seconds.count()
              << "s\n";
}

void fem_mesh::save_internal_variables(bool const have_converged)
{
    for (auto& submesh : submeshes)
    {
        submesh.save_internal_variables(have_converged);
    }
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
        std::string const& boundary_name = boundary["Name"];
        std::string const& boundary_type = boundary["Type"];

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
                           [&, dof_offset = std::ref(dof_offset)](auto const dof) {
                               return dof + dof_offset;
                           });

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
        for (auto const& boundary_variant : nonfollower_load.natural_interface())
        {
            std::visit(
                [&](auto const& boundary_mesh) {
                    auto const& times = boundary_mesh.time_history();

                    history.insert(begin(times), end(times));
                },
                boundary_variant);
        }
    }
    return {begin(history), end(history)};
}

void fem_mesh::check_boundary_conditions(json const& boundary_data) const
{
    for (auto const& boundary : boundary_data)
    {
        for (auto const& mandatory_field : {"Name", "Type"})
        {
            if (!boundary.count(mandatory_field))
            {
                throw std::domain_error("\"" + std::string(mandatory_field)
                                        + "\" was not specified in \"BoundaryCondition\".");
            }
        }
        if ((!boundary.count("Time")) && (!boundary.count("GenerateType")))
        {
            throw std::domain_error("Neither \"Time\" nor \"GenerateType\" was specified in "
                                    "\"BoundaryCondition\".");
        }
    }
}

void fem_mesh::write(std::int32_t const time_step, double const current_time)
{
    // nodal variables
    if (writer->is_output_requested("displacement"))
    {
        writer->field("displacement", coordinates->displacement(), 3);
    }
    if (writer->is_output_requested("reaction_force"))
    {
        writer->field("reaction_force", reaction_forces, 2);
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
                                                            coordinates->size() * 4,
                                                            convert(output),
                                                            output),
                                  4);
                }
            },
            output_variable);
    }
    writer->write(time_step, current_time);
}

void fem_mesh::allocate_variable_names()
{
    for (auto const& name : writer->outputs())
    {
        if (is_nodal_variable(name)) continue;

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
