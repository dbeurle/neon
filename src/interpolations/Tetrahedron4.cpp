
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
    // using NodalCoordinate = std::tuple<int, double, double, double>;
    //
    // // Initialize nodal coordinates array as r and s
    // std::array<NodalCoordinate, 4> constexpr local_coordinates{
    //     {{0, 1.0, 0.0, 0.0}, {1, 0.0, 1.0, 0.0}, {2, 0.0, 0.0, 1.0}, {3, 0.0, 0.0, 0.0}}};

    numerical_quadrature->evaluate([&](auto const& coordinate) {
        auto const& [l, r, s, t] = coordinate;

        vector N(4);
        matrix rhea(4, 3);

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

        return std::make_tuple(N, rhea);
    });

    extrapolation = matrix::Ones(nodes(), 1);
}
}
