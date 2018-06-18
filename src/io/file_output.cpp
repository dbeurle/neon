
#include "file_output.hpp"

#include <exception>

#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Winconsistent-missing-override"
#endif

#include "vtkCellData.h"
#include "vtkCellTypes.h"
#include "vtkDataObject.h"
#include "vtkInformation.h"

#if defined __clang__
#pragma clang diagnostic pop
#endif

#include <boost/filesystem.hpp>

namespace neon
{
namespace io
{
file_output::file_output(std::string file_name, json const& visualisation_data)
    : file_name(file_name), unstructured_mesh(vtkSmartPointer<vtkUnstructuredGrid>::New())
{
    if (visualisation_data.count("WriteEvery"))
    {
        write_every = visualisation_data["WriteEvery"];
    }

    boost::filesystem::path directory_path(directory_name);
    boost::filesystem::create_directory(directory_path);

    pvd_file.open(file_name + ".pvd");

    if (!pvd_file.is_open())
    {
        throw std::domain_error("Not able to write to disk for visualisation\n");
    }

    pvd_file << "<?xml version=\"1.0\"?>\n";
    pvd_file << "<VTKFile type=\"Collection\" version=\"0.1\">\n";
    pvd_file << std::string(2, ' ') << "<Collection>\n";

    unstructured_mesh->Allocate();

    if (visualisation_data.count("Fields"))
    {
        for (std::string const& field : visualisation_data["Fields"])
        {
            output_variables.insert(field);
        }
    }
}

file_output::~file_output()
{
    // Close off the last of the file for the timestepping
    pvd_file << std::string(2, ' ') << "</Collection>\n"
             << "</VTKFile>\n";
    pvd_file.close();
}

void file_output::write_to_file(int const time_step, double const total_time)
{
    auto unstructured_mesh_writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();

    auto const vtk_filename = directory_name + "/" + file_name + "_" + std::to_string(time_step)
                              + "." + unstructured_mesh_writer->GetDefaultFileExtension();

    std::cout << "\n"
              << std::string(4, ' ') << "Writing solution to file for step " << time_step << "\n";

    unstructured_mesh_writer->SetFileName(vtk_filename.c_str());
    unstructured_mesh_writer->SetInputData(unstructured_mesh);

    if (!use_binary_format) unstructured_mesh_writer->SetDataModeToAscii();

    unstructured_mesh_writer->Write();

    // Update the pvd file for timestep mapping
    pvd_file << std::string(4, ' ') << "<DataSet timestep = \"" << std::to_string(total_time)
             << "\" file = \"" << directory_name << "/" << file_name << "_"
             << std::to_string(time_step) << "."
             << unstructured_mesh_writer->GetDefaultFileExtension() << "\" />\n";
}

void file_output::add_field(std::string const& name, vector const& field, int const components)
{
    auto scalar_field = vtkSmartPointer<vtkDoubleArray>::New();

    scalar_field->SetName(name.c_str());
    scalar_field->SetNumberOfComponents(components);
    scalar_field->Allocate(field.size() / components);

    for (std::int64_t index{0}; index < field.size(); index += components)
    {
        scalar_field->InsertNextTuple(field.data() + index);
    }
    unstructured_mesh->GetPointData()->AddArray(scalar_field);
}

vtk_file_output::vtk_file_output(std::string file_name, json const& visualisation_data)
    : file_output(file_name, visualisation_data)
{
    pvd_file.open(file_name + ".pvd");

    if (!pvd_file.is_open())
    {
        throw std::domain_error("Not able to write to disk for visualisation\n");
    }

    pvd_file << "<?xml version=\"1.0\"?>\n";
    pvd_file << "<VTKFile type=\"Collection\" version=\"0.1\">\n";
    pvd_file << std::string(2, ' ') << "<Collection>\n";

    unstructured_mesh->Allocate();
}

vtk_file_output::~vtk_file_output()
{
    // close off the last of the file for the time stepping
    pvd_file << std::string(2, ' ') << "</Collection>\n"
             << "</VTKFile>\n";
    pvd_file.close();
}

void vtk_file_output::write(int const time_step, double const total_time)
{
    // wait on previous future
    write_future.wait();

    if (write_future.get() == 1)
    {
        throw std::domain_error("Error in VTK file IO occurred");
    }

    write_future = std::async(std::launch::async, [this, time_step, total_time]() {
        auto unstructured_mesh_writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();

        auto const vtk_filename = directory_name + "/" + file_name + "_" + std::to_string(time_step)
                                  + "." + unstructured_mesh_writer->GetDefaultFileExtension();

        std::cout << "\n"
                  << std::string(4, ' ') << "Writing solution to file for step " << time_step
                  << "\n";

        unstructured_mesh_writer->SetFileName(vtk_filename.c_str());
        unstructured_mesh_writer->SetInputData(unstructured_mesh);

        if (!use_binary_format) unstructured_mesh_writer->SetDataModeToAscii();

        // Update the pvd file for timestep mapping
        pvd_file << std::string(4, ' ') << "<DataSet timestep = \"" << std::to_string(total_time)
                 << "\" file = \"" << directory_name << "/" << file_name << "_"
                 << std::to_string(time_step) << "."
                 << unstructured_mesh_writer->GetDefaultFileExtension() << "\" />\n";

        return unstructured_mesh_writer->Write();
    });
}

void vtk_file_output::coordinates(matrix3x const& configuration) {}

void vtk_file_output::mesh(indices const& all_node_indices, element_topology const topology)
{
    auto const vtk_node_indices = convert_to_vtk(all_node_indices, topology);

    for (std::int64_t element{0}; element < vtk_node_indices.cols(); ++element)
    {
        auto node_indices = vtkSmartPointer<vtkIdList>::New();

        for (std::int64_t node{0}; node < vtk_node_indices.rows(); ++node)
        {
            node_indices->InsertNextId(static_cast<std::int64_t>(vtk_node_indices(node, element)));
        }
        unstructured_mesh->InsertNextCell(to_vtk(topology), node_indices);
    }
}

void vtk_file_output::field(std::string const& name,
                            vector const& flattened_field,
                            std::int64_t const components)
{
    // ensure not adding to a field still being written to file
    write_future.wait();

    auto vtk_field = vtkSmartPointer<vtkDoubleArray>::New();

    vtk_field->SetName(name.c_str());
    vtk_field->SetNumberOfComponents(components);
    vtk_field->Allocate(flattened_field.size() / components);

    for (std::int64_t index{0}; index < flattened_field.size(); index += components)
    {
        vtk_field->InsertNextTuple(flattened_field.data() + index);
    }
    unstructured_mesh->GetPointData()->AddArray(vtk_field);
}
}

namespace diffusion
{
file_output::file_output(std::string file_name, json const& visualisation_data, fem_mesh const& mesh)
    : io::file_output(file_name, visualisation_data), mesh(mesh)
{
    // Check the output set against the known values for this module
    for (auto const& output : output_variables)
    {
        if (vector_map.find(output) == vector_map.end() && output != primary_field)
        {
            throw std::domain_error("Output \"" + output
                                    + "\" is not valid for a heat diffusion simulation\n");
        }
    }
    add_mesh();
}

void file_output::write(int const time_step, double const total_time, vector const& scalars)
{
    if (time_step % write_every != 0) return;

    // Write out the required fields
    for (auto const& name : output_variables)
    {
        if (name == primary_field)
        {
            add_field(primary_field, scalars, 1);
        }
    }
    write_to_file(time_step, total_time);
}

void file_output::add_mesh()
{
    // Populate an unstructured grid object
    unstructured_mesh->SetPoints(io::vtk_coordinates(mesh.geometry().coordinates()));

    for (auto const& submesh : mesh.meshes())
    {
        auto const vtk_node_indices = convert_to_vtk(submesh.all_node_indices(), submesh.topology());

        for (std::int64_t element{0}; element < vtk_node_indices.cols(); ++element)
        {
            auto node_indices = vtkSmartPointer<vtkIdList>::New();

            for (std::int64_t node{0}; node < vtk_node_indices.rows(); ++node)
            {
                node_indices->InsertNextId(static_cast<std::int64_t>(vtk_node_indices(node, element)));
            }
            unstructured_mesh->InsertNextCell(to_vtk(submesh.topology()), node_indices);
        }
    }
}
}
}
