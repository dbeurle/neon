
#include "triangle.hpp"

#include "geometry/projection.hpp"

#include <Eigen/Geometry>

namespace neon
{
triangle3::triangle3(triangle_quadrature::point const p)
    : surface_interpolation(std::make_unique<triangle_quadrature>(p))
{
    this->precompute_shape_functions();
}

void triangle3::precompute_shape_functions()
{
    // Initialize nodal coordinates array as r and s
    m_quadrature->evaluate([&](auto const& coordinates) {
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

double triangle3::compute_measure(matrix const& nodal_coordinates) const
{
    // Use the cross product identity 2A = | a x b | to compute face area
    vector3 const direction0 = nodal_coordinates.col(0) - nodal_coordinates.col(2);
    vector3 const direction1 = nodal_coordinates.col(1) - nodal_coordinates.col(2);

    vector3 const normal = direction0.cross(direction1);

    return normal.norm() / 2.0;
}

triangle6::triangle6(triangle_quadrature::point const p)
    : surface_interpolation(std::make_unique<surface_quadrature>(triangle_quadrature(p)))
{
    this->precompute_shape_functions();
}

void triangle6::precompute_shape_functions()
{
    using NodalCoordinate = std::tuple<int, double, double>;

    // Initialize nodal coordinates array as Xi, Eta, Zeta
    std::array<NodalCoordinate, 6> constexpr local_coordinates{{
        {0, 1.0, 0.0},
        {1, 0.0, 1.0},
        {2, 0.0, 0.0},
        {3, 0.5, 0.5},
        {4, 0.0, 0.5},
        {5, 0.5, 0.0},
    }};

    matrix N_matrix(m_quadrature->points(), nodes());
    matrix local_quadrature_coordinates = matrix::Ones(m_quadrature->points(), 3);

    m_quadrature->evaluate([&](auto const& coordinate) {
        auto const& [l, r, s] = coordinate;

        auto const t = 1.0 - r - s;

        vector N(6);
        matrix rhea(6, 2);

        N(0) = r * (2.0 * r - 1.0);
        N(1) = s * (2.0 * s - 1.0);
        N(2) = t * (2.0 * t - 1.0);
        N(3) = 4.0 * r * s;
        N(4) = 4.0 * s * t;
        N(5) = 4.0 * r * t;

        // r coordinates
        rhea(0, 0) = 4.0 * r - 1.0;
        rhea(1, 0) = 0.0;
        rhea(2, 0) = -4.0 * t + 1.0;
        rhea(3, 0) = 4.0 * s;
        rhea(4, 0) = -4.0 * s;
        rhea(5, 0) = 4.0 * t - 4.0 * r;

        // s coordinates
        rhea(0, 1) = 0.0;
        rhea(1, 1) = 4.0 * s - 1.0;
        rhea(2, 1) = -4.0 * t + 1.0;
        rhea(3, 1) = 4.0 * r;
        rhea(4, 1) = 4.0 * t - 4.0 * s;
        rhea(5, 1) = -4.0 * r;

        local_quadrature_coordinates(l, 0) = r;
        local_quadrature_coordinates(l, 1) = s;

        N_matrix.row(l) = N;

        return std::make_tuple(N, rhea);
    });

    // Compute extrapolation algorithm matkrices
    matrix local_nodal_coordinates = matrix::Ones(nodes(), 3);

    for (auto const& [a, r, s] : local_coordinates)
    {
        local_nodal_coordinates(a, 0) = r;
        local_nodal_coordinates(a, 1) = s;
    }
    compute_extrapolation_matrix(N_matrix, local_nodal_coordinates, local_quadrature_coordinates);
}

double triangle6::compute_measure(matrix const& nodal_coordinates)
{
    return m_quadrature->integrate(0.0, [&](auto const& femval, auto const& l) {
        auto const& [N, dN] = femval;

        matrix2 const Jacobian = geometry::project_to_plane(nodal_coordinates) * dN;

        return Jacobian.determinant();
    });
}
}
