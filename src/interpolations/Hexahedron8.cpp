/*
 * neon - A finite element solver.
 *
 * For licensing please refer to the LICENSE.md file
 *
 */

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
        return std::make_tuple(N, rhea);
    });
}
}
