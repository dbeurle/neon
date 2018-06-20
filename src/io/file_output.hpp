
#pragma once

#include "numeric/index_types.hpp"
#include "numeric/dense_matrix.hpp"
#include "mesh/element_topology.hpp"

#include "io/json_forward.hpp"

#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>

#include <future>
#include <set>
#include <string>

namespace neon::io
{
class file_output
{
public:
    explicit file_output(std::string file_name, json const& visualisation_data);

    /// Virtual destructor to finish writing out the time step mapping
    virtual ~file_output() = default;

    /// Add coordinates to the file output set
    virtual void coordinates(matrix const& configuration) = 0;

    /// Add mesh information to the file output set
    virtual void mesh(indices const& all_node_indices, element_topology const topology) = 0;

    /// Write out to file in the format specified
    virtual void write(int const time_step, double const total_time) = 0;

    /// Add the field to the output field with
    /// \param name Field name
    /// \param data Flat vector of encoded data
    /// \param components Number of components encoded in the field
    virtual void field(std::string const& name, vector const& data, std::int64_t const components) = 0;

    [[nodiscard]] bool is_output_requested(std::string const& name) const;

    [[nodiscard]] std::set<std::string> outputs() const noexcept { return output_variables; }

protected:
    /// Directory to store visualisation output
    std::string const directory_name{"visualisation"};
    /// Output file name that appears on disk
    std::string file_name;
    /// Requested variables from the input file
    std::set<std::string> output_variables;
    /// Time steps to write out (e.g. two is every second time step)
    int write_every{1};
};

// vtk_file_output output should be asynchronous as these operations are
// typically expensive.  For the asynchronous aspect - the file_output class
// should take the vtk object to write out while the computation continues.
class vtk_file_output : public io::file_output
{
public:
    explicit vtk_file_output(std::string file_name, json const& visualisation_data);

    virtual ~vtk_file_output();

    /// Write out to file
    virtual void write(int const time_step, double const total_time) override final;

    virtual void coordinates(matrix const& configuration) override final;

    /// Add mesh information to the file output set
    virtual void mesh(indices const& all_node_indices, element_topology const topology) override final;

    virtual void field(std::string const& name,
                       vector const& data,
                       std::int64_t const components) override final;

private:
    /// Holds one if successful and zero otherwise.  Use this to perform the
    /// io without blocking computation
    std::future<int> future;
    /// VTK representation of mesh
    vtkSmartPointer<vtkUnstructuredGrid> unstructured_mesh;
    /// Default to using binary VTK output for efficiency
    bool use_binary_format{true};
    /// Stream for writing time history
    std::ofstream pvd_file;
};
}

// template <class fem_mesh>
// void file_output<fem_mesh>::write(int const time_step, double const total_time)
// {
//     if (time_step % write_every != 0) return;
//
//     vector nodal_averaged_value, insertions;
//
//     // Write out the required fields
//     for (auto const& name : output_variables)
//     {
//         if (auto found = tensor_map.find(name); found != tensor_map.end())
//         {
//             auto const& [name_str, name_enum] = *found;
//
//             auto const tensor_size = fem_mesh::internal_variable_type::tensor_size;
//
//             nodal_averaged_value.resize(mesh.geometry().size() * tensor_size);
//             insertions.resize(mesh.geometry().size() * tensor_size);
//
//             nodal_averaged_value.setZero();
//             insertions.setZero();
//
//             // Add internal variables
//             for (auto const& submesh : mesh.meshes())
//             {
//                 if (!submesh.internal_variables().has(name_enum))
//                 {
//                     throw std::domain_error("Internal variable " + name + " does not exist in mesh");
//                 }
//
//                 auto const [value, count] = submesh.nodal_averaged_variable(name_enum);
//
//                 nodal_averaged_value += value;
//
//                 insertions += count;
//             }
//             nodal_averaged_value = nodal_averaged_value.cwiseQuotient(insertions);
//
//             add_field(name_str, nodal_averaged_value, tensor_size);
//         }
//         else if (auto found = scalar_map.find(name); found != scalar_map.end())
//         {
//             auto const& [name_str, name_enum] = *found;
//
//             nodal_averaged_value.resize(mesh.geometry().size());
//             insertions.resize(mesh.geometry().size());
//
//             nodal_averaged_value.setZero();
//             insertions.setZero();
//
//             for (auto const& submesh : mesh.meshes())
//             {
//                 auto const [value, count] = submesh.nodal_averaged_variable(name_enum);
//                 nodal_averaged_value += value;
//                 insertions += count;
//             }
//             nodal_averaged_value = nodal_averaged_value.cwiseQuotient(insertions);
//
//             add_field(name_str, nodal_averaged_value, 1);
//         }
//         else
//         {
//             throw std::domain_error("Field \"" + name
//                                     + "\" was not found in mesh internal variables\n");
//         }
//     }
// }
