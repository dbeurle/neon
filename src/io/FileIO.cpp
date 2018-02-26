
#include "FileIO.hpp"

#include <exception>

#if defined __clang__
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
FileIO::FileIO(std::string file_name, json const& visualisation_data)
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
        for (auto const& field : visualisation_data["Fields"])
        {
            output_set.insert(field.get<std::string>());
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

void FileIO::add_field(std::string const& name, vector const& field, int const components)
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

namespace diffusion
{
FileIO::FileIO(std::string file_name, json const& visualisation_data, fem_mesh const& mesh)
    : neon::FileIO(file_name, visualisation_data), mesh(mesh)
{
    // Check the output set against the known values for this module
    for (auto const& output : output_set)
    {
        if (vector_map.find(output) == vector_map.end() && output != primary_field)
        {
            throw std::domain_error("Output \"" + output
                                    + "\" is not valid for a heat diffusion simulation\n");
        }
    }
    add_mesh();
}

void FileIO::write(int const time_step, double const total_time, vector const& scalars)
{
    if (time_step % write_every != 0) return;

    // Write out the required fields
    for (auto const& name : output_set)
    {
        if (name == primary_field)
        {
            add_field(primary_field, scalars, 1);
        }
    }
    write_to_file(time_step, total_time);
}

void FileIO::add_mesh()
{
    // Populate an unstructured grid object
    unstructured_mesh->SetPoints(io::vtk_coordinates(mesh.geometry().coordinates()));

    for (auto const& submesh : mesh.meshes())
    {
        auto const vtk_ordered_connectivity = convert_to_vtk(submesh.element_connectivity(),
                                                             submesh.topology());

        for (std::int64_t element{0}; element < vtk_ordered_connectivity.cols(); ++element)
        {
            auto vtk_node_list = vtkSmartPointer<vtkIdList>::New();

            for (std::int64_t node{0}; node < vtk_ordered_connectivity.rows(); ++node)
            {
                vtk_node_list->InsertNextId(
                    static_cast<std::int64_t>(vtk_ordered_connectivity(node, element)));
            }
            unstructured_mesh->InsertNextCell(to_vtk(submesh.topology()), vtk_node_list);
        }
    }
}
}
}
