
#pragma once

#include "constitutive/InternalVariables.hpp"
#include "mesh/NodeOrderingAdapter.hpp"
#include "mesh/diffusion/femMesh.hpp"
#include "mesh/solid/femMesh.hpp"

#include <map>
#include <string>
#include <unordered_set>

#include "vtkSmartPointer.h"

#include <json/forwards.h>

class vtkUnstructuredGrid;
class vtkXMLUnstructuredGridWriter;

namespace neon
{
class FileIO
{
public:
    explicit FileIO(std::string file_name, Json::Value const& visualisation_data);

    ~FileIO();

    FileIO(FileIO&&) = default;

protected:
    virtual void write_to_file(int const time_step, double const total_time);

    /** Write out the field to a vtk file */
    void add_field(std::string const& name, Vector const& data, int const components);

protected:
    std::string file_name;

    vtkSmartPointer<vtkUnstructuredGrid> unstructured_mesh;
    // vtkSmartPointer<vtkXMLUnstructuredGridWriter> unstructured_mesh_writer;

    std::ofstream pvd_file; //!< Stream for writing time history

    NodeOrderingAdapter adapter; //!< Converting neon to vtk ordering

    std::unordered_set<std::string> output_set;

    int write_every = 1; //!< Time steps to write out (e.g. two is every second time step)
    bool use_binary_format = false;
};

namespace solid
{
class FileIO : public neon::FileIO
{
public:
    using ScalarMap = std::map<std::string, InternalVariables::Scalar>;
    using TensorMap = std::map<std::string, InternalVariables::Tensor>;

    std::string const primary_field{"Displacement"};

public:
    explicit FileIO(std::string file_name,
                    Json::Value const& visualisation_data,
                    femMesh const& fem_mesh);

    FileIO(FileIO&&) = default;

    void write(int const time_step, double const total_time);

private:
    void add_mesh();

private:
    femMesh const& fem_mesh;

    // clang-format off
    ScalarMap const scalar_map{{"AccumulatedPlasticStrain", InternalVariables::Scalar::EffectivePlasticStrain},
                               {"VonMisesStress", InternalVariables::Scalar::VonMisesStress}};

    TensorMap const tensor_map{{"CauchyStress", InternalVariables::Tensor::Cauchy},
                               {"LinearisedStrain", InternalVariables::Tensor::LinearisedStrain},
                               {"LinearisedPlasticStrain", InternalVariables::Tensor::LinearisedPlasticStrain},
                               {"DeformationGradient", InternalVariables::Tensor::DeformationGradient},
                               {"DisplacementGradient", InternalVariables::Tensor::DisplacementGradient}};
    // clang-format on
};
}

namespace diffusion
{
class FileIO : public neon::FileIO
{
public:
    using VectorMap = std::map<std::string, InternalVariables::Vector>;

    std::string const primary_field{"Temperature"};

public:
    explicit FileIO(std::string file_name,
                    Json::Value const& visualisation_data,
                    femMesh const& fem_mesh);

    ~FileIO();

    FileIO(FileIO&&) = default;

    void write(int const time_step, double const total_time, Vector const& temperature);

private:
    void add_mesh();

private:
    femMesh const& fem_mesh;

    VectorMap const vector_map{{"HeatFlux", InternalVariables::Vector::HeatFlux}};
};
}
}
