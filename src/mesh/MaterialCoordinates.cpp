
#include "MaterialCoordinates.hpp"

#include <range/v3/algorithm/for_each.hpp>

#include <range/v3/view/iota.hpp>
#include <range/v3/view/zip.hpp>

#include "vtkDoubleArray.h"
#include "vtkPoints.h"

namespace neon
{
template <int coordinate, int stride>
using vector_view = Eigen::Map<vector, coordinate, Eigen::InnerStride<stride>>;

MaterialCoordinates::MaterialCoordinates(Vector const& initial_coordinates)
    : NodalCoordinates(initial_coordinates), x(initial_coordinates)
{
}

Vector MaterialCoordinates::displacement(List const& local_dofs) const
{
    using namespace ranges;

    Vector localdisp(local_dofs.size());

    for_each(view::zip(view::ints(0), local_dofs), [&](auto const& zip_pair) {
        auto const & [i, local_dof] = zip_pair;
        localdisp(i) = x(local_dof) - X(local_dof);
    });
    return localdisp;
}

void MaterialCoordinates::update_current_xy_configuration(vector const& u)
{
    std::cout << "Hallo.....!" << std::endl;

    std::cout << "Size of the vector " << u.size() << std::endl;
    std::cout << "Number of nodes " << u.size() / 2 << std::endl;

    std::cout << "Coordinates of the vector " << x.size() << std::endl;
    std::cout << "Number of nodes " << x.size() / 3 << std::endl;

    vector_view<0, 3>(x.data(),
                      x.size() / 3) = vector_view<0, 3>(X.data(), X.size() / 3)
                                      + vector_view<0, 2>(const_cast<vector::Scalar*>(u.data()),
                                                          u.size() / 2);

    vector_view<1, 3>(x.data(),
                      x.size() / 3) = vector_view<1, 3>(X.data(), X.size() / 3)
                                      + vector_view<1, 2>(const_cast<vector::Scalar*>(u.data()),
                                                          u.size() / 2);
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
    auto displacements = vtkSmartPointer<vtkDoubleArray>::New();
    displacements->SetNumberOfComponents(3);
    displacements->SetName("Displacements");

    for (auto i = 0; i < X.size(); i += 3)
    {
        displacements->InsertNextTuple3(x(i) - X(i), x(i + 1) - X(i + 1), x(i + 2) - X(i + 2));
    }
    return displacements;
}

matrix3x MaterialCoordinates::get_configuration(local_indices const& local_nodes,
                                                vector const& configuration) const
{
    matrix3x element_displacement(3, local_nodes.size());

    for (auto lnode = 0; lnode < local_nodes.size(); lnode++)
    {
        element_displacement.col(lnode) = configuration.segment<3>(3 * local_nodes[lnode]);
    }
    return element_displacement;
}
}
