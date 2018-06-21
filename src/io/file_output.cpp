
#include "file_output.hpp"

#include "mesh/node_ordering_adapter.hpp"
#include "io/json.hpp"

#include <exception>

#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Winconsistent-missing-override"
#endif

#include "io/vtk_coordinates.hpp"

#include <vtkIdList.h>
#include <vtkDoubleArray.h>
#include <vtkIdTypeArray.h>
#include <vtkPointData.h>
#include <vtkXMLUnstructuredGridWriter.h>

#if defined __clang__
#pragma clang diagnostic pop
#endif

#include <boost/filesystem.hpp>

namespace neon::io
{
file_output::file_output(std::string file_name, json const& visualisation_data)
    : file_name(file_name)
{
    if (visualisation_data.is_null())
    {
        throw std::domain_error("Visualisation data must be specified");
    }
    if (file_name.empty())
    {
        throw std::domain_error("Name field must be specified");
    }

    if (visualisation_data.find("WriteEvery") != visualisation_data.end())
    {
        write_every = visualisation_data["WriteEvery"];
    }
    if (visualisation_data.find("Fields") != visualisation_data.end())
    {
        for (std::string const& field : visualisation_data["Fields"])
        {
            output_variables.insert(field);
        }
    }
    else
    {
        std::cout << std::string(4, ' ') << "No outputs were requested.  I find this strange.\n";
    }
    boost::filesystem::create_directory(boost::filesystem::path(directory_name));
}

bool file_output::is_output_requested(std::string const& name) const
{
    return output_variables.find(name) != end(output_variables);
}

vtk_file_output::vtk_file_output(std::string file_name, json const& visualisation_data)
    : file_output(file_name, visualisation_data),
      unstructured_mesh(vtkSmartPointer<vtkUnstructuredGrid>::New())
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
    future.wait();

    // close off the last of the file for the time stepping
    pvd_file << std::string(2, ' ') << "</Collection>\n"
             << "</VTKFile>\n";
    pvd_file.close();
}

void vtk_file_output::write(int const time_step, double const current_time)
{
    // wait on previous future
    if (future.valid())
    {
        if (future.wait(); future.get() == 0)
        {
            throw std::domain_error("Error in VTK file IO occurred");
        }
    }

    future = std::async(std::launch::async, [this, time_step, current_time]() {
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
        pvd_file << std::string(4, ' ') << "<DataSet timestep = \"" << std::to_string(current_time)
                 << "\" file = \"" << directory_name << "/" << file_name << "_"
                 << std::to_string(time_step) << "."
                 << unstructured_mesh_writer->GetDefaultFileExtension() << "\" />\n";

        return unstructured_mesh_writer->Write();
    });
}

void vtk_file_output::coordinates(matrix const& configuration)
{
    auto points = vtkSmartPointer<vtkPoints>::New();

    points->Allocate(configuration.size() / configuration.rows());

    for (std::int64_t i{0}; i < configuration.cols(); ++i)
    {
        points->InsertNextPoint(configuration(0, i),
                                (configuration.rows() > 1 ? configuration(1, i) : 0.0),
                                (configuration.rows() > 2 ? configuration(2, i) : 0.0));
    }
    unstructured_mesh->SetPoints(points);
}

void vtk_file_output::mesh(indices const& all_node_indices, element_topology const topology)
{
    indices const vtk_node_indices = convert_to_vtk(all_node_indices, topology);

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
                            vector const& field_vector,
                            std::int64_t const components)
{
    if (future.valid())
    {
        future.wait();
    }

    auto vtk_field = vtkSmartPointer<vtkDoubleArray>::New();

    vtk_field->SetName(name.c_str());
    vtk_field->SetNumberOfComponents(components);
    vtk_field->Allocate(field_vector.size() / components);

    for (std::int64_t index{0}; index < field_vector.size(); index += components)
    {
        vtk_field->InsertNextTuple(field_vector.data() + index);
    }
    unstructured_mesh->GetPointData()->AddArray(vtk_field);
}
}
