
#include "mesh.hpp"

#include "math/block_sequence.hpp"
#include "mesh/basic_mesh.hpp"
#include "mesh/dof_allocator.hpp"
#include "io/post/variable_string_adapter.hpp"
#include "io/json.hpp"

#include <chrono>
#include <exception>
#include <memory>
#include <numeric>

#include <termcolor/termcolor.hpp>

namespace neon::mechanics::beam
{
static bool is_nodal_variable(std::string const& name)
{
    return name == "displacement" || name == "rotation";
}

mesh::mesh(basic_mesh const& basic_mesh,
           json const& material_data,
           json const& simulation_data,
           double const generate_time_step,
           std::map<std::string, std::unique_ptr<geometry::profile>> const& profile_store)
    : coordinates(std::make_shared<material_coordinates>(basic_mesh.coordinates())),
      displacement(active_dofs() / 2),
      rotation(active_dofs() / 2),
      generate_time_step{generate_time_step},
      writer(std::make_unique<io::vtk_file_output>(simulation_data["Name"],
                                                   simulation_data["Visualisation"]))
{
    check_boundary_conditions(simulation_data["BoundaryConditions"]);

    writer->coordinates(coordinates->coordinates());

    std::cout << "simulation name is: " << simulation_data["Name"] << std::endl;

    for (auto const& section : simulation_data["sections"])
    {
        if (section.find("name") == section.end())
        {
            throw std::domain_error("A section is missing a \"name\" field");
        }

        for (auto const& submesh : basic_mesh.meshes(section["name"].get<std::string>()))
        {
            submeshes.emplace_back(material_data, simulation_data, section, coordinates, submesh);

            writer->mesh(submesh.all_node_indices(), submesh.topology());
        }
    }

    allocate_boundary_conditions(simulation_data, basic_mesh);
    allocate_variable_names();
}

bool mesh::is_symmetric() const { return true; }

void mesh::update_internal_variables(vector const& displacement_rotation, double const time_step_size)
{
    auto const start = std::chrono::steady_clock::now();

    displacement = displacement_rotation(block_sequence<3, 6>{0, displacement.size()});
    rotation = displacement_rotation(block_sequence<3, 6>{3, rotation.size()});

    coordinates->update_current_configuration(displacement);

    for (auto& submesh : submeshes)
    {
        submesh.update_internal_variables(time_step_size);
    }

    auto const end = std::chrono::steady_clock::now();
    std::chrono::duration<double> const elapsed_seconds = end - start;

    std::cout << std::string(6, ' ') << "Internal variable update took " << elapsed_seconds.count()
              << "s\n";
}

bool mesh::is_nonfollower_load(std::string const& boundary_type) const
{
    return boundary_type == "moment" || boundary_type == "force";
}

void mesh::allocate_boundary_conditions(json const& simulation_data, basic_mesh const& basic_mesh)
{
    // Populate the boundary conditions and their corresponding mesh
    for (auto const& boundary_data : simulation_data["BoundaryConditions"])
    {
        std::string const& boundary_name = boundary_data["Name"];
        std::string const& boundary_type = boundary_data["Type"];

        if (boundary_type == "displacement" || boundary_type == "rotation")
        {
            this->allocate_dirichlet_boundary(boundary_type, boundary_data, basic_mesh);
        }
        else if (is_nonfollower_load(boundary_type))
        {
            nonfollower_loads.emplace(boundary_name,
                                      nonfollower_load_boundary(coordinates,
                                                                basic_mesh.meshes(boundary_name),
                                                                simulation_data,
                                                                boundary_data,
                                                                dof_table,
                                                                generate_time_step));
        }
        else
        {
            throw std::domain_error("BoundaryCondition \"" + boundary_type + "\" is not recognised");
        }
    }
}

void mesh::allocate_dirichlet_boundary(std::string const& boundary_type,
                                       json const& boundary,
                                       basic_mesh const& basic_mesh)
{
    std::string const& boundary_name = boundary["Name"];

    for (auto const& [dof_key, dof_offset] : dof_table)
    {
        if (boundary.find(dof_key) != boundary.end())
        {
            auto boundary_dofs = unique_dof_allocator<traits::dofs_per_node>(
                basic_mesh.meshes(boundary_name));

            // Offset the degrees of freedom on the boundary
            std::transform(begin(boundary_dofs),
                           end(boundary_dofs),
                           begin(boundary_dofs),
                           [&, dof_offset = std::ref(dof_offset)](auto const dof) {
                               return dof + dof_offset + (boundary_type == "displacement" ? 0 : 3);
                           });

            m_dirichlet_boundaries[boundary_name].emplace_back(boundary_dofs,
                                                               boundary,
                                                               dof_key,
                                                               generate_time_step);
        }
    }
}

std::vector<double> mesh::time_history() const
{
    std::set<double> history;

    // Append time history from each boundary condition
    for (auto const& [key, boundaries] : m_dirichlet_boundaries)
    {
        for (auto const& boundary : boundaries)
        {
            auto const times = boundary.time_history();
            history.insert(begin(times), end(times));
        }
    }
    // for (auto const& [key, nonfollower_load] : nonfollower_loads)
    // {
    //     for (auto const& boundary_variant : nonfollower_load.natural_interface())
    //     {
    //         std::visit(
    //             [&](auto const& boundary_mesh) {
    //                 auto const times = boundary_mesh.time_history();
    //                 history.insert(begin(times), end(times));
    //             },
    //             boundary_variant);
    //     }
    // }
    return {begin(history), end(history)};
}

void mesh::write(std::int32_t const time_step, double const current_time)
{
    // nodal variables
    if (writer->is_output_requested("displacement"))
    {
        writer->field("displacement", displacement, 3);
    }
    if (writer->is_output_requested("rotation"))
    {
        writer->field("rotation", rotation, 3);
    }
    writer->write(time_step, current_time);
}

void mesh::check_boundary_conditions(json const& boundary_data) const
{
    for (auto const& boundary : boundary_data)
    {
        for (auto const& mandatory_field : {"Name", "Type"})
        {
            if (boundary.find(mandatory_field) == boundary.end())
            {
                throw std::domain_error("\"" + std::string(mandatory_field)
                                        + "\" was not specified in \"BoundaryCondition\".");
            }
        }
        if (boundary.find("Time") == boundary.end() && boundary.find("GenerateType") == boundary.end())
        {
            throw std::domain_error("Neither \"Time\" nor \"GenerateType\" was specified in "
                                    "\"BoundaryCondition\".");
        }
    }
}

void mesh::allocate_variable_names()
{
    for (auto const& name : writer->outputs())
    {
        if (is_nodal_variable(name)) continue;

        output_variables.emplace_back(variable::convert(name));
    }
}
}
