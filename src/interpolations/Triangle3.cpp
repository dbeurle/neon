
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
    // using NodalCoordinate = std::tuple<int, double, double>;

    // Initialize nodal coordinates array as r and s
    // std::array<NodalCoordinate, 3> constexpr local_coordinates = {{
    //     {0, 1.0, 0.0},
    //     {1, 0.0, 1.0},
    //     {2, 0.0, 0.0},
    // }};

    numerical_quadrature->evaluate([&](auto const& coordinates) {
        auto const& [l, r, s] = coordinates;

        vector N(3);
        matrix rhea(3, 2);

        N(0) = r;
        N(1) = s;
        N(2) = 1.0 - r - s;

        rhea(0, 0) = 1.0;
        rhea(1, 0) = 0.0;
        rhea(2, 0) = -1.0;

        rhea(0, 1) = 0.0;
        rhea(1, 1) = 1.0;
        rhea(2, 1) = -1.0;

        return std::make_tuple(N, rhea);
    });

    extrapolation = matrix::Ones(nodes(), 1);
}

double Triangle3::compute_measure(matrix const& nodal_coordinates) const
{
    // Use the cross product identity 2A = | a x b | to compute face area
    vector3 const direction0 = nodal_coordinates.col(0) - nodal_coordinates.col(2);
    vector3 const direction1 = nodal_coordinates.col(1) - nodal_coordinates.col(2);

    vector3 normal = direction0.cross(direction1);

    return normal.norm() / 2.0;
}
}
