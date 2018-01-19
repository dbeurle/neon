
#include "MaterialCoordinates.hpp"

#include "vtkDoubleArray.h"
#include "vtkPoints.h"

namespace neon
{
MaterialCoordinates::MaterialCoordinates(matrix3x const& initial_coordinates)
    : NodalCoordinates(initial_coordinates), x(initial_coordinates)
{
}

matrix3x MaterialCoordinates::displacement() const { return x - X; }

void MaterialCoordinates::update_current_xy_configuration(vector const& u)
{
    x.row(0) = X.row(0) + u.transpose()(Eigen::seq(0, u.size() - 1, 2));
    x.row(1) = X.row(1) + u.transpose()(Eigen::seq(1, u.size() - 1, 2));
}

void MaterialCoordinates::update_current_configuration(vector const& u)
{
    x.row(0) = X.row(0) + u.transpose()(Eigen::seq(0, u.size() - 1, 3));
    x.row(1) = X.row(1) + u.transpose()(Eigen::seq(1, u.size() - 1, 3));
    x.row(2) = X.row(2) + u.transpose()(Eigen::seq(2, u.size() - 1, 3));
}

vtkSmartPointer<vtkPoints> MaterialCoordinates::vtk_coordinates() const
{
    auto points = vtkSmartPointer<vtkPoints>::New();

    points->Allocate(X.size() / 3);

    for (auto i = 0; i < X.size(); i += 3)
    {
        points->InsertNextPoint(X(i), X(i + 1), X(i + 2));
    }
    return points;
}

vtkSmartPointer<vtkDoubleArray> MaterialCoordinates::vtk_displacement() const
{
    static_assert(!X.IsRowMajor, "This assumes the storage is column major");

    auto displacements = vtkSmartPointer<vtkDoubleArray>::New();
    displacements->Allocate(X.cols());
    displacements->SetNumberOfComponents(X.rows());
    displacements->SetName("Displacements");

    for (auto i = 0; i < X.cols(); i++)
    {
        displacements->InsertNextTuple3(x(0, i) - X(0, i), x(1, i) - X(1, i), x(2, i) - X(2, i));
    }
    return displacements;
}
}
