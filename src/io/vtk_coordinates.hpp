
#include "vtkDoubleArray.h"
#include "vtkPoints.h"
#include "vtkSmartPointer.h"

#pragma once

namespace neon::io
{
/**
 * Copy coordinates from a matrix into a vtkPoints data structure
 * @param X matrix of coordinates
 * @return A vktSmartPointer of vtkPoints from \p X
 */
template <typename MatrixType>
vtkSmartPointer<vtkPoints> vtk_coordinates(MatrixType const& X)
{
    auto points = vtkSmartPointer<vtkPoints>::New();

    points->Allocate(X.size() / X.rows());

    for (auto i = 0; i < X.size(); i += X.rows())
    {
        points->InsertNextPoint(X(i),
                                (X.rows() > 1 ? X(i + 1) : 0.0),
                                (X.rows() > 2 ? X(i + 2) : 0.0));
    }
    return points;
}

/**
 * Copy displacements into a vtkDoubleArray for output
 * @return vtkSmartPointer containing a vtkDoubleArray
 */
template <typename MatrixType>
vtkSmartPointer<vtkDoubleArray> vtk_displacement(MatrixType const& u)
{
    // Enforce column major storage
    static_assert(!u.IsRowMajor, "This assumes the storage is column major");

    // The coordinates should be known at compile time for performance reasons
    static_assert(MatrixType::RowsAtCompileTime != -1, "The number of rows must greater than one.");
    // static_assert(u.RowsAtCompileTime > 3, "The number of rows must be less than three");

    auto displacements = vtkSmartPointer<vtkDoubleArray>::New();
    displacements->Allocate(u.cols());
    displacements->SetNumberOfComponents(u.rows());
    displacements->SetName("Displacements");

    std::array<double, u.RowsAtCompileTime> displacement_tuple;

    for (auto i = 0; i < u.cols(); i++)
    {
        for (auto d = 0; d < u.rows(); ++d)
        {
            displacement_tuple[d] = u(d, i);
        }
        displacements->InsertNextTuple(displacement_tuple.data());
    }
    return displacements;
}
}
