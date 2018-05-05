
#include "vtkDoubleArray.h"
#include "vtkPoints.h"
#include "vtkSmartPointer.h"

#pragma once

namespace neon::io
{
/// Copy coordinates from a matrix into a vtkPoints data structure
/// \param X matrix of coordinates
/// \return A vktSmartPointer of vtkPoints from \p X
template <typename matrix_type>
vtkSmartPointer<vtkPoints> vtk_coordinates(matrix_type const& X)
{
    auto points = vtkSmartPointer<vtkPoints>::New();

    points->Allocate(X.size() / X.rows());

    for (std::int64_t i{0}; i < X.size(); i += X.rows())
    {
        points->InsertNextPoint(X(i),
                                (X.rows() > 1 ? X(i + 1) : 0.0),
                                (X.rows() > 2 ? X(i + 2) : 0.0));
    }
    return points;
}

/// Copy input data into a vtkDoubleArray for output
/// \tparam rows Compile time rows of input data
/// \param input_data Data to output size = rows * cols
/// \param array_name Name of the resulting array
/// \return vtkSmartPointer containing a vtkDoubleArray
template <int rows>
vtkSmartPointer<vtkDoubleArray> vtk_displacement(matrixdx<rows> const& input_data,
                                                 std::string const array_name = "displacement")
{
    static_assert(rows != Eigen::Dynamic, "The number of rows must be compile-time fixed");

    static_assert(!input_data.IsRowMajor, "This assumes the storage is column major");

    static_assert(std::is_same<typename matrixxd<rows>::value_type, double>::value,
                  "Input data must be doubles");

    auto vtk_double_array = vtkSmartPointer<vtkDoubleArray>::New();

    vtk_double_array->Allocate(input_data.cols());
    vtk_double_array->SetNumberOfComponents(rows);
    vtk_double_array->SetName(array_name.c_str());

    for (std::int64_t i{0}; i < input_data.cols(); i++)
    {
        vtk_double_array->InsertNextTuple(input_data.col(i).data());
    }
    return vtk_double_array;
}

/// Copy input data into a vtkDoubleArray for output where the input data
/// is encoded
/// \param nodal_variable Nodal variable to use
/// \param start_index Variable to start from
/// \param stride Size of the components
/// \param array_name Pretty name to appear in ParaView
/// \return vtkSmartPointer containing a vtkDoubleArray
template <typename nodal_variable_type>
vtkSmartPointer<vtkDoubleArray> vtk_array(nodal_variable_type const& nodal_variable,
                                          std::int64_t const start_index,
                                          std::int64_t const stride,
                                          std::string const array_name)
{
    static_assert(std::is_same<typename nodal_variable_type::value_type, double>::value,
                  "Input data must be doubles");

    auto vtk_double_array = vtkSmartPointer<vtkDoubleArray>::New();

    vtk_double_array->Allocate(nodal_variable.size());
    vtk_double_array->SetNumberOfComponents(stride);
    vtk_double_array->SetName(array_name.c_str());

    for (std::int64_t i{0}; i < nodal_variable.size() * nodal_variable.components; i += stride)
    {
        vtk_double_array->InsertNextTuple(nodal_variable.data().col(i).data() + start_index);
    }
    return vtk_double_array;
}
}
