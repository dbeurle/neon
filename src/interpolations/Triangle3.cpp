
#include "Triangle3.hpp"

#include <Eigen/Geometry>

namespace neon
{
Triangle3::Triangle3(TriangleQuadrature::Rule rule)
    : SurfaceInterpolation(std::make_unique<TriangleQuadrature>(rule))
{
    this->precompute_shape_functions();
}

void Triangle3::precompute_shape_functions()
{
    using NodalCoordinate = std::tuple<int, double, double>;

    // Initialize nodal coordinates array as Xi, Eta, Zeta
    std::array<NodalCoordinate, 3> constexpr local_coordinates = {{
        {0, 0.0, 0.0},
        {1, 1.0, 0.0},
        {2, 0.0, 1.0},
    }};

    Matrix N_matrix(numerical_quadrature->points(), nodes());
    Matrix local_quadrature_coordinates = Matrix::Ones(numerical_quadrature->points(), 3);

    numerical_quadrature->evaluate([&](auto const& coordinates) {
        auto const& [l, r, s] = coordinates;

        Vector N(3);
        Matrix rhea(3, 2);

        N(0) = r;
        N(1) = s;
        N(2) = 1.0 - r - s;

        rhea(0, 0) = 1.0;
        rhea(1, 0) = 0.0;
        rhea(2, 0) = -1.0;

        rhea(0, 1) = 0.0;
        rhea(1, 1) = 1.0;
        rhea(2, 1) = -1.0;

        local_quadrature_coordinates(l, 0) = r;
        local_quadrature_coordinates(l, 1) = s;

        N_matrix.row(l) = N;

        return std::make_tuple(N, rhea);
    });

    // Compute extrapolation algorithm matrices
    Matrix local_nodal_coordinates = Matrix::Ones(nodes(), 3);

    for (auto const& [a, r, s] : local_coordinates)
    {
        local_nodal_coordinates(a, 0) = r;
        local_nodal_coordinates(a, 1) = s;
    }
    compute_extrapolation_matrix(N_matrix, local_nodal_coordinates, local_quadrature_coordinates);
}

double Triangle3::compute_measure(Matrix const& nodal_coordinates) const
{
    // Use the cross product identity 2A = | a x b | to compute face area
    Vector3 const direction0 = nodal_coordinates.col(0) - nodal_coordinates.col(2);
    Vector3 const direction1 = nodal_coordinates.col(1) - nodal_coordinates.col(2);

    Vector3 normal = direction0.cross(direction1);

    return normal.norm() / 2.0;
}
}
