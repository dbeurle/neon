
#include "Tetrahedron4.hpp"

namespace neon
{
Tetrahedron4::Tetrahedron4(TetrahedronQuadrature::Rule rule)
    : VolumeInterpolation(std::make_unique<TetrahedronQuadrature>(rule))
{
    this->precompute_shape_functions();
}

void Tetrahedron4::precompute_shape_functions()
{
    using NodalCoordinate = std::tuple<int, double, double, double>;

    // Initialize nodal coordinates array as Xi, Eta, Zeta
    constexpr std::array<NodalCoordinate, 4> local_coordinates = {{
        {0, 0.0, 0.0, 0.0},
        {1, 1.0, 0.0, 0.0},
        {2, 0.0, 1.0, 0.0},
        {3, 0.0, 0.0, 1.0},
    }};

    Matrix N_matrix(numerical_quadrature->points(), nodes());
    Matrix local_quadrature_coordinates = Matrix::Ones(numerical_quadrature->points(), 4);

    numerical_quadrature->evaluate([&](auto const& coordinate) {
        auto const & [ l, r, s, t ] = coordinate;

        Vector N(4);
        Matrix rhea(4, 3);

        N(0) = r;
        N(1) = s;
        N(2) = t;
        N(3) = 1.0 - r - s - t;

        rhea(0, 0) = 1.0;
        rhea(0, 1) = 0.0;
        rhea(0, 2) = 0.0;

        rhea(1, 0) = 0.0;
        rhea(1, 1) = 1.0;
        rhea(1, 2) = 0.0;

        rhea(2, 0) = 0.0;
        rhea(2, 1) = 0.0;
        rhea(2, 2) = 1.0;

        rhea(3, 0) = -1.0;
        rhea(3, 1) = -1.0;
        rhea(3, 2) = -1.0;

        local_quadrature_coordinates(l, 0) = r;
        local_quadrature_coordinates(l, 1) = s;
        local_quadrature_coordinates(l, 2) = t;

        N_matrix.row(l) = N;

        return std::make_tuple(N, rhea);
    });

    // Compute extrapolation algorithm matrices
    Matrix local_nodal_coordinates = Matrix::Ones(nodes(), 4);

    for (auto const & [ a, xi_a, eta_a, zeta_a ] : local_coordinates)
    {
        local_nodal_coordinates(a, 0) = xi_a;
        local_nodal_coordinates(a, 1) = eta_a;
        local_nodal_coordinates(a, 2) = zeta_a;
    }
    compute_extrapolation_matrix(N_matrix,
                                 local_nodal_coordinates,
                                 local_quadrature_coordinates);
}
}
