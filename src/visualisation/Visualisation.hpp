
#pragma once

#include "constitutive/InternalVariables.hpp"

#include "mesh/NodeOrderingAdapter.hpp"

#include <map>
#include <unordered_set>

#include "vtkSmartPointer.h"

#include <json/forwards.h>

class vtkUnstructuredGrid;
class vtkXMLUnstructuredGridWriter;

namespace neon
{
// Forward declaration of the continuum mechanics finite element mesh
namespace solid
{
class femMesh;
}

class Visualisation
{
public:
    Visualisation(std::string file_name,
                  solid::femMesh const& fem_mesh,
                  Json::Value const& visualisation_data);

    ~Visualisation();

    Visualisation(Visualisation&&) = default;

    void write(int time_step, double total_time);

protected:
    void allocate_field_maps();

    /** Fill geometry and initial result structures */
    void allocate_static_mesh();

    void write_tensor_field(std::string const& pretty_name,
                            InternalVariables::Tensor const& tensor_enum);

    void write_scalar_field(std::string const& pretty_name,
                            InternalVariables::Scalar const& scalar_enum);

    void write_primary_field();

protected:
    std::string file_name;

    solid::femMesh const& fem_mesh;

    vtkSmartPointer<vtkUnstructuredGrid> unstructured_mesh;
    // vtkSmartPointer<vtkXMLUnstructuredGridWriter> unstructured_mesh_writer;
    std::ofstream pvd_file;

    NodeOrderingAdapter adapter; //!< Converting neon to vtk ordering

    std::unordered_set<std::string> requested_fields;

    int write_every = 1; //!< Time steps to write out (e.g. two is every second time step)
    bool use_binary_format = true;

    std::map<std::string, InternalVariables::Tensor> string_to_tensor;
    std::map<std::string, InternalVariables::Scalar> string_to_scalar;
};
}
