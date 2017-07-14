
#include "mesh/solid/femMesh.hpp"

#include "mesh/BasicMesh.hpp"

#include <exception>

#include <memory>
#include <numeric>

#include <json/json.h>
#include <termcolor/termcolor.hpp>

#include <range/v3/action.hpp>
#include <range/v3/view.hpp>

#include "vtkCellData.h"
#include "vtkCellTypes.h"
#include "vtkDataObject.h"
#include "vtkDoubleArray.h"
#include "vtkIdList.h"
#include "vtkIdTypeArray.h"
#include "vtkInformation.h"
#include "vtkInformationQuadratureSchemeDefinitionVectorKey.h"
#include "vtkPointData.h"
#include "vtkPolyData.h"
#include "vtkPolyDataMapper.h"
#include "vtkQuadraturePointsGenerator.h"
#include "vtkQuadratureSchemeDefinition.h"
#include "vtkUnstructuredGrid.h"
#include "vtkUnstructuredGridReader.h"
#include "vtkXMLUnstructuredGridReader.h"
#include "vtkXMLUnstructuredGridWriter.h"

#include "vtkSmartPointer.h"

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

    allocate_boundary_conditions(simulation_data["BoundaryConditions"], basic_mesh);
}

femMesh::~femMesh() { finalise_vtk(); }

int femMesh::active_dofs() const { return 3 * material_coordinates->size(); }

void femMesh::internal_restart(Json::Value const& simulation_data)
{
    if (simulation_data["BoundaryConditions"].empty())
    {
        for (auto & [ name, boundaries ] : dirichlet_boundaries)
        {
            std::cout << termcolor::yellow << std::string(2, ' ') << "Boundary conditions for \""
                      << name << "\" have been inherited from the last load step"
                      << termcolor::reset << std::endl;

            for (auto& boundary : boundaries) boundary.inherit_from_last();
        }
        return;
    }

    auto const& boundary_data = simulation_data["BoundaryConditions"];

    check_boundary_conditions(boundary_data);
    reallocate_boundary_conditions(boundary_data);
}

void femMesh::update_internal_variables(Vector const& u, double const Δt)
{
    material_coordinates->update_current_configuration(u);

    for (auto& submesh : submeshes) submesh.update_internal_variables(Δt);
}

void femMesh::save_internal_variables(bool const have_converged)
{
    for (auto& submesh : submeshes) submesh.save_internal_variables(have_converged);
}

void femMesh::write(int const time_step, double const time)
{
    // Create an unstructured grid object
    auto unstructured_grid = vtkSmartPointer<vtkUnstructuredGrid>::New();
    unstructured_grid->Allocate();

    last_time_step = time_step;
    time_history.push_back(time);

    unstructured_grid->SetPoints(material_coordinates->vtk_coordinates());

    for (auto const& submesh : submeshes)
    {
        for (auto const& node_list : submesh.connectivities())
        {
            auto vtk_node_list = vtkSmartPointer<vtkIdList>::New();
            for (auto const& node : node_list)
            {
                vtk_node_list->InsertNextId(static_cast<long>(node));
            }
            unstructured_grid->InsertNextCell(adapter.to_vtk(submesh.topology()), vtk_node_list);
        }
    }

    unstructured_grid->GetPointData()->AddArray(material_coordinates->vtk_displacement());

    {
        Vector nodal_averaged_value = Vector::Zero(material_coordinates->size() * 9);
        Vector running_count = Vector::Zero(material_coordinates->size() * 9);

        // Add internal variables
        for (auto const& submesh : submeshes)
        {
            auto const[value, count] =
                submesh.nodal_averaged_variable(InternalVariables::Tensor::Cauchy);

            nodal_averaged_value += value;
            running_count += count;
        }

        // Average nodal values
        nodal_averaged_value = nodal_averaged_value.cwiseQuotient(running_count);

        // Put this into a vtkDoubleArray
        auto tensor_value = vtkSmartPointer<vtkDoubleArray>::New();

        tensor_value->SetNumberOfComponents(9);
        tensor_value->SetName("Cauchy stress");

        for (auto i = 0; i < material_coordinates->size() * 9; i += 9)
        {
            tensor_value->InsertNextTuple9(nodal_averaged_value(i + 0),
                                           nodal_averaged_value(i + 1),
                                           nodal_averaged_value(i + 2),
                                           nodal_averaged_value(i + 3),
                                           nodal_averaged_value(i + 4),
                                           nodal_averaged_value(i + 5),
                                           nodal_averaged_value(i + 6),
                                           nodal_averaged_value(i + 7),
                                           nodal_averaged_value(i + 8));
        }
        unstructured_grid->GetPointData()->AddArray(tensor_value);
    }

    auto unstructured_grid_writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();

    auto const filename = "Test_" + std::to_string(time_step) + ".vtu";

    unstructured_grid_writer->SetFileName(filename.c_str());
    unstructured_grid_writer->SetInputData(unstructured_grid);
    unstructured_grid_writer->SetDataModeToAscii();
    unstructured_grid_writer->Write();
}

void femMesh::finalise_vtk() const
{
    using namespace ranges;

    std::ofstream pvd_file;
    pvd_file.open("Test.pvd");

    pvd_file << "<?xml version=\"1.0\"?>\n";
    pvd_file << "<VTKFile type=\"Collection\" version=\"0.1\">\n";
    pvd_file << std::string(2, ' ') << "<Collection>\n";

    for_each(view::zip(time_history, view::ints(0)), [&](auto const& tpl) {
        auto const & [ time, i ] = tpl;
        pvd_file << std::string(4, ' ') << "<DataSet timestep = \"" << std::to_string(time)
                 << "\" file = \"Test_" << std::to_string(i) << ".vtu\" />\n";
    });

    pvd_file << std::string(2, ' ') << "</Collection>\n";
    pvd_file << "</VTKFile>\n";

    pvd_file.close();
}

void femMesh::allocate_boundary_conditions(Json::Value const& boundary_data,
                                           BasicMesh const& basic_mesh)
{
    using namespace ranges;

    // Populate the boundary meshes
    for (auto const& boundary : boundary_data)
    {
        auto const& boundary_name = boundary["Name"].asString();

        if (boundary["Type"].asString() == "Displacement")
        {
            auto const dirichlet_dofs = filter_dof_list(basic_mesh.meshes(boundary_name));

            for (auto const& name : boundary["Values"].getMemberNames())
            {
                auto const& dof_offset = dof_table.find(name)->second;

                dirichlet_boundaries[boundary_name]
                    .emplace_back(view::transform(dirichlet_dofs,
                                                  [&](auto dof) { return dof + dof_offset; }),
                                  boundary["Values"][name].asDouble());
            }
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
                for (auto& dirichlet_boundary : dirichlet_boundaries[boundary_name])
                {
                    dirichlet_boundary.update_value(boundary["Values"][name].asDouble());
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
                           [](auto const& submesh) { return submesh.connectivities(); }) |
           action::join | action::join | action::sort | action::unique |
           action::transform([=](auto const& i) { return i * 3; });
}
}
