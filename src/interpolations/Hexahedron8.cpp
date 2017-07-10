
#include "Hexahedron8.hpp"

#include <array>
#include <memory>
#include <tuple>

namespace neon
{
Hexahedron8::Hexahedron8(HexahedronQuadrature::Rule rule)
    : VolumeInterpolation(std::make_unique<HexahedronQuadrature>(rule))
{
    this->precompute_shape_functions();
}

void Hexahedron8::precompute_shape_functions()
{
    using NodalCoordinate = std::tuple<int, double, double, double>;

    // Initialize nodal coordinates array as Xi, Eta, Zeta
    constexpr std::array<NodalCoordinate, 8> local_coordinates = {{{0, -1.0, -1.0, -1.0},
                                                                   {1, 1.0, -1.0, -1.0},
                                                                   {2, 1.0, 1.0, -1.0},
                                                                   {3, -1.0, 1.0, -1.0},
                                                                   {4, -1.0, -1.0, 1.0},
                                                                   {5, 1.0, -1.0, 1.0},
                                                                   {6, 1.0, 1.0, 1.0},
                                                                   {7, -1.0, 1.0, 1.0}}};

    Matrix N_matrix(numerical_quadrature->points(), nodes());
    Matrix local_quadrature_coordinates = Matrix::Ones(numerical_quadrature->points(), 4);

    numerical_quadrature->evaluate([&](auto const& coordinate) {
        auto const & [ l, xi, eta, zeta ] = coordinate;

        Vector N(8);
        Matrix rhea(8, 3);

        for (auto const & [ a, xi_a, eta_a, zeta_a ] : local_coordinates)
        {
            N(a) = 1.0 / 8.0 * (1.0 + xi_a * xi) * (1.0 + eta_a * eta) * (1.0 + zeta_a * zeta);
            rhea(a, 0) = 1.0 / 8.0 * xi_a * (1.0 + eta_a * eta) * (1.0 + zeta_a * zeta);
            rhea(a, 1) = 1.0 / 8.0 * (1.0 + xi_a * xi) * eta_a * (1.0 + zeta_a * zeta);
            rhea(a, 2) = 1.0 / 8.0 * (1.0 + xi_a * xi) * (1.0 + eta_a * eta) * zeta_a;
        }
        local_quadrature_coordinates(l, 0) = xi;
        local_quadrature_coordinates(l, 1) = eta;
        local_quadrature_coordinates(l, 2) = zeta;

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
                                 local_quadrature_coordinates,
                                 nodes());
}
}
