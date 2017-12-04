
#include "FileIO.hpp"

#include "mesh/mech/solid/femMesh.hpp"

#include <json/value.h>

#include <exception>

#pragma clang diagnostic ignored "-Winconsistent-missing-override"

#include "vtkCellData.h"
#include "vtkCellTypes.h"
#include "vtkDataObject.h"
#include "vtkDoubleArray.h"
#include "vtkIdList.h"
#include "vtkIdTypeArray.h"
#include "vtkInformation.h"
#include "vtkInformationQuadratureSchemeDefinitionVectorKey.h"
#include "vtkPointData.h"
#include "vtkQuadraturePointsGenerator.h"
#include "vtkQuadratureSchemeDefinition.h"
#include "vtkUnstructuredGrid.h"
#include "vtkUnstructuredGridReader.h"
#include "vtkXMLUnstructuredGridReader.h"
#include "vtkXMLUnstructuredGridWriter.h"

#include <boost/filesystem.hpp>

namespace neon
{
FileIO::FileIO(std::string file_name, Json::Value const& visualisation_data)
    : file_name(file_name), unstructured_mesh(vtkSmartPointer<vtkUnstructuredGrid>::New())
{
    if (visualisation_data.isMember("WriteEvery"))
    {
        write_every = visualisation_data["WriteEvery"].asInt();
    }

    boost::filesystem::path directory_path(directory_name);
    boost::filesystem::create_directory(directory_path);

    pvd_file.open(file_name + ".pvd");

    if (!pvd_file.is_open())
    {
        throw std::runtime_error("Not able to write to disk for visualisation\n");
    }

    pvd_file << "<?xml version=\"1.0\"?>\n";
    pvd_file << "<VTKFile type=\"Collection\" version=\"0.1\">\n";
    pvd_file << std::string(2, ' ') << "<Collection>\n";

    unstructured_mesh->Allocate();

    if (visualisation_data.isMember("Fields"))
    {
        for (auto const& field : visualisation_data["Fields"])
        {
            output_set.insert(field.asString());
        }
    }
}

FileIO::~FileIO()
{
    // Close off the last of the file for the timestepping
    pvd_file << std::string(2, ' ') << "</Collection>\n"
             << "</VTKFile>\n";
    pvd_file.close();
}

void FileIO::write_to_file(int const time_step, double const total_time)
{
    auto unstructured_mesh_writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();

    auto const vtk_filename = directory_name + "/" + file_name + "_" + std::to_string(time_step)
                              + "." + unstructured_mesh_writer->GetDefaultFileExtension();

    std::cout << "\n"
              << std::string(4, ' ') << "Writing solution to file for step " << time_step << "\n";

    unstructured_mesh_writer->SetFileName(vtk_filename.c_str());
    unstructured_mesh_writer->SetInputData(unstructured_mesh);

    if (use_binary_format == false) unstructured_mesh_writer->SetDataModeToAscii();

    unstructured_mesh_writer->Write();

    // Update the pvd file for timestep mapping
    pvd_file << std::string(4, ' ') << "<DataSet timestep = \"" << std::to_string(total_time)
             << "\" file = \"" << directory_name << "/" << file_name << "_"
             << std::to_string(time_step) << "."
             << unstructured_mesh_writer->GetDefaultFileExtension() << "\" />\n";
}

void FileIO::add_field(std::string const& name, Vector const& field, int const components)
{
    auto scalar_value = vtkSmartPointer<vtkDoubleArray>::New();

    scalar_value->SetName(name.c_str());
    scalar_value->SetNumberOfComponents(components);
    // scalar_value->SetNumberOfTuples(field.size() / components);

    for (auto i = 0; i < field.size(); i += components)
    {
        scalar_value->InsertNextTuple(&field(i));
    }
    unstructured_mesh->GetPointData()->AddArray(scalar_value);
}

namespace mech::solid
{
FileIO::FileIO(std::string file_name, Json::Value const& visualisation_data, femMesh const& fem_mesh)
    : neon::FileIO(file_name, visualisation_data), fem_mesh(fem_mesh)
{
    // Check the output set against the known values for this module
    for (auto const& output : output_set)
    {
        if (scalar_map.find(output) == scalar_map.end()
            && tensor_map.find(output) == tensor_map.end() && output != primary_field)
        {
            throw std::runtime_error("Output \"" + output
                                     + "\" is not valid for a solid mechanics simulation\n");
        }
    }
    add_mesh();
}

void FileIO::write(int const time_step, double const total_time)
{
    if (time_step % write_every != 0) return;

    Vector nodal_averaged_value, running_count;

    // Write out the required fields
    for (auto const& name : output_set)
    {
        if (auto found = tensor_map.find(name); found != tensor_map.end())
        {
            nodal_averaged_value.resize(fem_mesh.coordinates().size() * 9);
            running_count.resize(fem_mesh.coordinates().size() * 9);

            nodal_averaged_value.setZero();
            running_count.setZero();

            // Add internal variables
            for (auto const& submesh : fem_mesh.meshes())
            {
                if (!submesh.internal_variables().has(found->second))
                {
                    throw std::runtime_error("Internal variable " + name + " does not exist in mesh");
                }
                auto const [value, count] = submesh.nodal_averaged_variable(found->second);
                nodal_averaged_value += value;
                running_count += count;
            }
            nodal_averaged_value = nodal_averaged_value.cwiseQuotient(running_count);

            add_field(found->first, nodal_averaged_value, 9);
        }
        else if (auto found = scalar_map.find(name); found != scalar_map.end())
        {
            nodal_averaged_value.resize(fem_mesh.coordinates().size());
            running_count.resize(fem_mesh.coordinates().size());

            nodal_averaged_value.setZero();
            running_count.setZero();

            for (auto const& submesh : fem_mesh.meshes())
            {
                auto const [value, count] = submesh.nodal_averaged_variable(found->second);
                nodal_averaged_value += value;
                running_count += count;
            }
            nodal_averaged_value = nodal_averaged_value.cwiseQuotient(running_count);

            add_field(found->first, nodal_averaged_value, 1);
        }
        else if (name == primary_field)
        {
            unstructured_mesh->GetPointData()->AddArray(fem_mesh.coordinates().vtk_displacement());
        }
        else
        {
            throw std::runtime_error("Field \"" + name
                                     + "\" was not found in mesh internal variables\n");
        }
    }
    write_to_file(time_step, total_time);
}

void FileIO::add_mesh()
{
    // Populate an unstructured grid object
    unstructured_mesh->SetPoints(fem_mesh.coordinates().vtk_coordinates());

    for (auto const& submesh : fem_mesh.meshes())
    {
        auto const vtk_ordered_connectivity = adapter.convert_to_vtk(submesh.connectivities(),
                                                                     submesh.topology());
        for (auto const& node_list : vtk_ordered_connectivity)
        {
            auto vtk_node_list = vtkSmartPointer<vtkIdList>::New();

            for (auto const& node : node_list)
            {
                vtk_node_list->InsertNextId(static_cast<long>(node));
            }
            unstructured_mesh->InsertNextCell(adapter.to_vtk(submesh.topology()), vtk_node_list);
        }
    }
    unstructured_mesh->GetPointData()->AddArray(fem_mesh.coordinates().vtk_displacement());
}
}

namespace diffusion
{
FileIO::FileIO(std::string file_name, Json::Value const& visualisation_data, femMesh const& fem_mesh)
    : neon::FileIO(file_name, visualisation_data), fem_mesh(fem_mesh)
{
    // Check the output set against the known values for this module
    for (auto const& output : output_set)
    {
        if (vector_map.find(output) == vector_map.end() && output != primary_field)
        {
            throw std::runtime_error("Output \"" + output
                                     + "\" is not valid for a heat diffusion simulation\n");
        }
    }
    add_mesh();
}

void FileIO::write(int const time_step, double const total_time, Vector const& scalars)
{
    if (time_step % write_every != 0) return;

    // Write out the required fields
    for (auto const& name : output_set)
    {
        if (name == primary_field) add_field(primary_field, scalars, 1);
    }
    write_to_file(time_step, total_time);
}

void FileIO::add_mesh()
{
    // Populate an unstructured grid object
    unstructured_mesh->SetPoints(fem_mesh.coordinates().vtk_coordinates());

    for (auto const& submesh : fem_mesh.meshes())
    {
        auto const vtk_ordered_connectivity = adapter.convert_to_vtk(submesh.connectivities(),
                                                                     submesh.topology());
        for (auto const& node_list : vtk_ordered_connectivity)
        {
            auto vtk_node_list = vtkSmartPointer<vtkIdList>::New();

            for (auto const& node : node_list)
            {
                vtk_node_list->InsertNextId(static_cast<long>(node));
            }
            unstructured_mesh->InsertNextCell(adapter.to_vtk(submesh.topology()), vtk_node_list);
        }
    }
}
}
}
