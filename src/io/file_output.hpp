
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

// namespace neon::mechanical
// {
// template <typename fem_mesh>
// class file_output : public io::file_output
// {
// public:
//     /// Mechanical mesh type
//     using fem_mesh_type = fem_mesh;
//
//     using variable_type = typename fem_mesh_type::internal_variable_type;
//
//     using scalar_map_t = std::map<std::string, typename variable_type::scalar>;
//     using tensor_map_t = std::map<std::string, typename variable_type::second>;
//
// public:
//     std::string const primary_field{"Displacement"};
//
// public:
//     explicit file_output(std::string file_name,
//                          json const& visualisation_data,
//                          fem_mesh_type const& mesh);
//
//     file_output(file_output&&) = default;
//
//     void write(int const time_step, double const total_time);
//
// private:
//     void add_mesh();
//
// private:
//     /// Reference to the mesh to output
//     fem_mesh_type const& mesh;
//     /// Allowable fields must also be in either scalar_map / tensor_map
//     std::set<std::string> allowable_fields{{"AccumulatedPlasticStrain",
//                                             "VonMisesStress",
//                                             "Damage",
//                                             "EnergyReleaseRate",
//                                             "CauchyStress",
//                                             "LinearisedStrain",
//                                             "LinearisedPlasticStrain",
//                                             "DeformationGradient",
//                                             "DisplacementGradient",
//                                             "KinematicHardening",
//                                             "BackStress",
//                                             "Displacement",
//                                             "ReactionForce",
//                                             "ActiveShearModulus",
//                                             "InactiveShearModulus",
//                                             "ActiveSegments",
//                                             "InactiveSegments",
//                                             "ReductionFactor"}};
//
//     // clang-format off
//     scalar_map_t const scalar_map{{"AccumulatedPlasticStrain", variable_type::scalar::effective_plastic_strain},
//                                   {"VonMisesStress", variable_type::scalar::von_mises_stress},
//                                   {"Damage", variable_type::scalar::damage},
//                                   {"ActiveShearModulus", variable_type::scalar::active_shear_modulus},
//                                   {"InactiveShearModulus", variable_type::scalar::inactive_shear_modulus},
//                                   {"ActiveSegments", variable_type::scalar::active_segments},
//                                   {"InactiveSegments", variable_type::scalar::inactive_segments},
//                                   {"ReductionFactor", variable_type::scalar::reduction_factor},
//                                   {"EnergyReleaseRate", variable_type::scalar::energy_release_rate}};
//
//     tensor_map_t const tensor_map{{"CauchyStress", variable_type::second::cauchy_stress},
//                                   {"LinearisedStrain", variable_type::second::linearised_strain},
//                                   {"LinearisedPlasticStrain", variable_type::second::linearised_plastic_strain},
//                                   {"DeformationGradient", variable_type::second::deformation_gradient},
//                                   {"DisplacementGradient", variable_type::second::displacement_gradient},
//                                   {"KinematicHardening", variable_type::second::kinematic_hardening},
//                                   {"BackStress", variable_type::second::back_stress}};
//     // clang-format on
// };
//
// template <class fem_mesh>
// file_output<fem_mesh>::file_output(std::string file_name,
//                                    json const& visualisation_data,
//                                    fem_mesh const& mesh)
//     : io::file_output(file_name, visualisation_data), mesh(mesh)
// {
//     // Check the output set against the known values for this module
//     if (!std::includes(begin(allowable_fields),
//                        end(allowable_fields),
//                        begin(output_variables),
//                        end(output_variables)))
//     {
//         throw std::domain_error("Requested output is not valid for a solid mechanics "
//                                 "simulation\n");
//     }
//     add_mesh();
// }
//
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
//         else if (name == primary_field)
//         {
//             unstructured_mesh->GetPointData()->AddArray(
//                 io::vtk_displacement(mesh.geometry().displacement()));
//         }
//         else if (name == "ReactionForce")
//         {
//             add_field("reaction forces", mesh.nodal_reaction_forces(), fem_mesh::traits::size);
//         }
//         else
//         {
//             throw std::domain_error("Field \"" + name
//                                     + "\" was not found in mesh internal variables\n");
//         }
//     }
//     write_to_file(time_step, total_time);
// }
//
// template <class fem_mesh>
// void file_output<fem_mesh>::add_mesh()
// {
//     // Populate an unstructured grid object
//     unstructured_mesh->SetPoints(io::vtk_coordinates(mesh.geometry().coordinates()));
//
//     for (auto const& submesh : mesh.meshes())
//     {
//         auto const vtk_node_indices = convert_to_vtk(submesh.all_node_indices(), submesh.topology());
//
//         for (std::int64_t element{0}; element < vtk_node_indices.cols(); ++element)
//         {
//             auto node_indices = vtkSmartPointer<vtkIdList>::New();
//
//             for (std::int64_t node{0}; node < vtk_node_indices.rows(); ++node)
//             {
//                 node_indices->InsertNextId(static_cast<std::int64_t>(vtk_node_indices(node, element)));
//             }
//             unstructured_mesh->InsertNextCell(to_vtk(submesh.topology()), node_indices);
//         }
//     }
//     unstructured_mesh->GetPointData()->AddArray(io::vtk_displacement(mesh.geometry().coordinates()));
// }
// }
