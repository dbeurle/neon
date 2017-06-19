/*
 * neon - A finite element solver.
 *
 * For licensing please refer to the LICENSE.md file
 *
 */

#include "Tetrahedron10.hpp"

#include <memory>

namespace neon
{
Tetrahedron10::Tetrahedron10(TetrahedronQuadrature::Rule rule)
    : VolumeInterpolation(std::make_unique<VolumeQuadrature>(TetrahedronQuadrature(rule)))
{
    this->precompute_shape_functions();
}

void Tetrahedron10::precompute_shape_functions()
{
    numerical_quadrature->evaluate([&](auto const& coordinate) {
        auto const & [ l, r, s, t ] = coordinate;

        auto const u = 1.0 - r - s - t;

        Vector N(10);
        Matrix rhea(10, 3);

        N(0) = r * (2.0 * r - 1.0);
        N(1) = s * (2.0 * s - 1.0);
        N(2) = t * (2.0 * t - 1.0);
        N(3) = u * (2.0 * u - 1.0);
        N(4) = 4.0 * r * s;
        N(5) = 4.0 * s * t;
        N(6) = 4.0 * t * u;
        N(7) = 4.0 * r * u;
        N(8) = 4.0 * r * t;
        N(9) = 4.0 * s * u;

        rhea(0, 0) = 4.0 * r - 1.0;
        rhea(1, 0) = 0.0;
        rhea(2, 0) = 0.0;
        rhea(3, 0) = -3.0 + 4.0 * r + 4.0 * s + 4.0 * t;
        rhea(4, 0) = 4.0 * s;
        rhea(5, 0) = 0.0;
        rhea(6, 0) = -4.0 * t;
        rhea(7, 0) = 4.0 - 4.0 * s - 8.0 * r - 4.0 * t;
        rhea(8, 0) = 4.0 * t;
        rhea(9, 0) = -4.0 * s;

        rhea(0, 1) = 0.0;
        rhea(1, 1) = 4.0 * s - 1.0;
        rhea(2, 1) = 0.0;
        rhea(3, 1) = -3.0 + 4.0 * r + 4.0 * s + 4.0 * t;
        rhea(4, 1) = 4.0 * r;
        rhea(5, 1) = 4.0 * t;
        rhea(6, 1) = -4.0 * t;
        rhea(7, 1) = -4.0 * r;
        rhea(8, 1) = 0.0;
        rhea(9, 1) = 4.0 - 4.0 * r - 8.0 * s - 4.0 * t;

        rhea(0, 2) = 0.0;
        rhea(1, 2) = 0.0;
        rhea(2, 2) = 4.0 * t - 1.0;
        rhea(3, 2) = -3.0 + 4.0 * r + 4.0 * s + 4.0 * t;
        rhea(4, 2) = 0.0;
        rhea(5, 2) = 4.0 * s;
        rhea(6, 2) = 4.0 - 4.0 * r - 4.0 * s - 8.0 * t;
        rhea(7, 2) = -4.0 * r;
        rhea(8, 2) = 4.0 * r;
        rhea(9, 2) = -4.0 * s;

        return std::make_tuple(N, rhea);
    });
}
}
