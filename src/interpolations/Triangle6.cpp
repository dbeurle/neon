
#include "Triangle6.hpp"

namespace neon
{
Triangle6::Triangle6(TriangleQuadrature::Rule rule)
    : SurfaceInterpolation(std::make_unique<SurfaceQuadrature>(TriangleQuadrature(rule)))
{
    this->precompute_shape_functions();
}

void Triangle6::precompute_shape_functions()
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

    matrix N_matrix(numerical_quadrature->points(), nodes());
    matrix local_quadrature_coordinates = matrix::Ones(numerical_quadrature->points(), 3);

    numerical_quadrature->evaluate([&](auto const& coordinate) {
        auto const& [l, r, s] = coordinate;

        auto const t = 1.0 - r - s;

        Vector N(6);
        Matrix rhea(6, 2);

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

double Triangle6::compute_measure(Matrix const& nodal_coordinates)
{
    // Use numerical quadrature to compute int () dA
    double face_area = 0.0;

    // auto const L = numerical_quadrature->points();
    //
    // Matrix rhea(6, 2);
    // Matrix planarCoordinates = this->projectCoordinatesToPlane(nodal_coordinates);

    // for (int l = 0; l < L; l++)
    // {
    //     rhea.col(0) = dN_dXi.col(l);
    //     rhea.col(1) = dN_dEta.col(l);
    //
    //     Eigen::Matrix2d Jacobian = planarCoordinates * rhea;
    //
    //     face_area += Jacobian.determinant() * quadratureWeight(l);
    // }
    return face_area / 2.0;
}
}
