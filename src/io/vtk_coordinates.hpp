
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

    for (std::int64_t i = 0; i < X.size(); i += X.rows())
    {
        points->InsertNextPoint(X(i),
                                (X.rows() > 1 ? X(i + 1) : 0.0),
                                (X.rows() > 2 ? X(i + 2) : 0.0));
    }
    return points;
}

/// Copy input data into a vtkDoubleArray for output
/// \return vtkSmartPointer containing a vtkDoubleArray
template <typename matrix_type>
vtkSmartPointer<vtkDoubleArray> vtk_displacement(matrix_type const& input_data,
                                                 std::string const array_name = "displacement")
{
    static_assert(!input_data.IsRowMajor, "This assumes the storage is column major");

    static_assert(std::is_same<typename matrix_type::value_type, double>::value,
                  "Input data must be doubles");

    static_assert(matrix_type::RowsAtCompileTime != Eigen::Dynamic,
                  "The number of rows must be compile-time fixed");

    auto vtk_double_array = vtkSmartPointer<vtkDoubleArray>::New();

    vtk_double_array->Allocate(input_data.cols());
    vtk_double_array->SetNumberOfComponents(input_data.rows());
    vtk_double_array->SetName(array_name.c_str());

    for (std::int64_t i = 0; i < input_data.cols(); i++)
    {
        vtk_double_array->InsertNextTuple(input_data.col(i).data());
    }
    return vtk_double_array;
}
}
