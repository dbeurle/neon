
#include "mesh/solid/femMesh.hpp"

#include "mesh/BasicMesh.hpp"

#include <exception>
#include <json/json.h>
#include <memory>
#include <numeric>

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
    auto const& simulation_name = simulation_data["Name"].asString();

    // Populate the volume meshes
    for (auto const& submesh : basic_mesh.meshes(simulation_name))
    {
        submeshes.emplace_back(material_data, simulation_data, material_coordinates, submesh);
    }

    // Populate the boundary meshes
    for (auto const& boundary : simulation_data["BoundaryConditions"])
    {
        if (boundary["Name"].empty())
        {
            throw std::runtime_error("Missing Name in BoundaryConditions\n");
        }
        if (boundary["Type"].empty())
        {
            throw std::runtime_error("Missing Type in BoundaryConditions\n");
        }

        auto const& boundary_name = boundary["Name"].asString();
        auto const& boundary_type = boundary["Type"].asString();
        auto const& boundary_meshes = basic_mesh.meshes(boundary_name);

        if (boundary_type == "Displacement")
        {
            for (auto const& name : boundary["Values"].getMemberNames())
            {
                using namespace ranges;

                auto const dof_offset = dof_table.find(name)->second;

                // Filter the data by collecting all the connectivities,
                // placing these into a flat array, finding the unique entries
                // and finally offsetting for the correct nodal dof
                List const dirichlet_dofs =
                    view::transform(boundary_meshes,
                                    [](auto const& submesh) { return submesh.connectivities(); }) |
                    action::join | action::join | action::sort | action::unique |
                    action::transform([=](auto const& i) { return i * 3 + dof_offset; });

                dirichlet_boundaries[boundary_name].emplace_back(dirichlet_dofs,
                                                                 boundary["Values"][name].asDouble());
            }
        }
    }
}

int femMesh::active_dofs() const { return 3 * material_coordinates->size(); }

void femMesh::update_internal_variables(Vector const& u)
{
    material_coordinates->update_current_configuration(u);

    for (auto& submesh : submeshes) submesh.update_internal_variables();
}

void femMesh::save_internal_variables(bool const have_converged)
{
    for (auto& submesh : submeshes) submesh.save_internal_variables(have_converged);
}

void femMesh::write(int filename_append) const
{
    // Create an unstructured grid object
    auto unstructured_grid = vtkSmartPointer<vtkUnstructuredGrid>::New();
    unstructured_grid->Allocate();

    unstructured_grid->SetPoints(material_coordinates->vtk_coordinates());

    NodeOrderingAdapter adapter;

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

    auto const filename =
        filename_append < 0 ? "Test.vtu" : "Test_" + std::to_string(filename_append) + ".vtu";

    unstructured_grid_writer->SetFileName(filename.c_str());
    unstructured_grid_writer->SetInputData(unstructured_grid);
    unstructured_grid_writer->SetDataModeToAscii();
    unstructured_grid_writer->Write();
}
}
