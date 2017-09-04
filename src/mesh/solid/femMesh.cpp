
#include "mesh/solid/femMesh.hpp"

#include "mesh/BasicMesh.hpp"

#include <exception>
#include <memory>
#include <numeric>

#include <json/json.h>
#include <termcolor/termcolor.hpp>

#include <range/v3/action.hpp>
#include <range/v3/view.hpp>

namespace neon::solid
{
femMesh::femMesh(BasicMesh const& basic_mesh,
                 Json::Value const& material_data,
                 Json::Value const& simulation_data)
    : material_coordinates(std::make_shared<MaterialCoordinates>(basic_mesh.coordinates()))
{
    check_boundary_conditions(simulation_data["BoundaryConditions"]);

    auto const& simulation_name = simulation_data["Name"].asString();

    for (auto const& submesh : basic_mesh.meshes(simulation_name))
    {
        submeshes.emplace_back(material_data, simulation_data, material_coordinates, submesh);
    }

    allocate_boundary_conditions(simulation_data, basic_mesh);
}

void femMesh::internal_restart(Json::Value const& simulation_data)
{
    if (!simulation_data.isMember("BoundaryConditions"))
    {
        for (auto & [ name, boundaries ] : displacement_bcs)
        {
            std::cout << termcolor::yellow << std::string(2, ' ') << "Boundary conditions for \""
                      << name << "\" have been inherited from the last load step"
                      << termcolor::reset << std::endl;

            for (auto& boundary : boundaries) boundary.internal_restart();
        }
    }
    else
    {
        auto const& boundary_data = simulation_data["BoundaryConditions"];
        check_boundary_conditions(boundary_data);
        reallocate_boundary_conditions(boundary_data);
    }
}

void femMesh::update_internal_variables(Vector const& u, double const time_step_size)
{
    auto start = std::chrono::high_resolution_clock::now();

    material_coordinates->update_current_configuration(u);

    for (auto& submesh : submeshes) submesh.update_internal_variables(time_step_size);

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;

    std::cout << std::string(6, ' ') << "Internal variable update took " << elapsed_seconds.count()
              << "s\n";
}

void femMesh::save_internal_variables(bool const have_converged)
{
    for (auto& submesh : submeshes) submesh.save_internal_variables(have_converged);
}

void femMesh::allocate_boundary_conditions(Json::Value const& simulation_data,
                                           BasicMesh const& basic_mesh)
{
    using namespace ranges;

    auto const& boundary_data = simulation_data["BoundaryConditions"];

    // Populate the boundary meshes
    for (auto const& boundary : boundary_data)
    {
        auto const& boundary_name = boundary["Name"].asString();

        if (auto const& boundary_type = boundary["Type"].asString(); boundary_type == "Displacement")
        {
            auto const dirichlet_dofs = filter_dof_list(basic_mesh.meshes(boundary_name));

            for (auto const& name : boundary["Values"].getMemberNames())
            {
                auto const& dof_offset = dof_table.find(name)->second;

                if (dof_table.find(name) == dof_table.end())
                {
                    throw std::runtime_error("x, y or z are acceptable "
                                             "coordinates\n");
                }

                displacement_bcs[boundary_name].emplace_back(view::transform(dirichlet_dofs,
                                                                             [&](auto const& dof) {
                                                                                 return dof
                                                                                        + dof_offset;
                                                                             }),
                                                             boundary["Values"][name].asDouble());
            }
        }
        else if (boundary_type == "Traction")
        {
            for (auto const& name : boundary["Values"].getMemberNames())
            {
                if (dof_table.find(name) == dof_table.end())
                {
                    throw std::runtime_error("x, y or z are acceptable "
                                             "coordinates\n");
                }

                auto const& dof_offset = dof_table.find(name)->second;

                nf_loads[boundary_name].emplace_back(material_coordinates,
                                                     basic_mesh.meshes(boundary_name),
                                                     dof_offset,
                                                     boundary["Values"][name].asDouble(),
                                                     simulation_data);
            }
        }
        else
        {
            throw std::runtime_error("BoundaryCondition \"" + boundary_type + "\" is not recognised");
        }
    }
}

void femMesh::reallocate_boundary_conditions(Json::Value const& boundary_data)
{
    using namespace ranges;

    for (auto const& boundary : boundary_data)
    {
        auto const& boundary_name = boundary["Name"].asString();

        if (boundary["Type"].asString() == "Displacement")
        {
            for (auto const& name : boundary["Values"].getMemberNames())
            {
                for (auto& dirichlet_boundary : displacement_bcs[boundary_name])
                {
                    dirichlet_boundary.internal_restart(boundary["Values"][name].asDouble());
                }
            }
        }
    }
}

void femMesh::check_boundary_conditions(Json::Value const& boundary_data) const
{
    for (auto const& boundary : boundary_data)
    {
        if (boundary["Name"].empty())
        {
            throw std::runtime_error("Missing \"Name\" in BoundaryConditions\n");
        }
        if (boundary["Type"].empty())
        {
            throw std::runtime_error("Missing \"Type\" in BoundaryConditions\n");
        }
    }
}

List femMesh::filter_dof_list(std::vector<SubMesh> const& boundary_mesh) const
{
    using namespace ranges;

    return view::transform(boundary_mesh,
                           [](auto const& submesh) { return submesh.connectivities(); })
           | action::join | action::join | action::sort | action::unique
           | action::transform([=](auto const& i) { return i * 3; });
}
}
