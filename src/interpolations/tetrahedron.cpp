
#include "tetrahedron.hpp"

namespace neon
{
tetrahedron4::tetrahedron4(tetrahedron_quadrature::Rule rule)
    : volume_interpolation(std::make_unique<tetrahedron_quadrature>(rule))
{
    this->precompute_shape_functions();
}

void tetrahedron4::precompute_shape_functions()
{
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

tetrahedron10::tetrahedron10(tetrahedron_quadrature::Rule rule)
    : volume_interpolation(std::make_unique<volume_quadrature>(tetrahedron_quadrature(rule)))
{
    this->precompute_shape_functions();
}

void tetrahedron10::precompute_shape_functions()
{
    using NodalCoordinate = std::tuple<int, double, double, double>;

    constexpr std::array<NodalCoordinate, 10> local_coordinates = {{{0, 1.0, 0.0, 0.0},
                                                                    {1, 0.0, 1.0, 0.0},
                                                                    {2, 0.0, 0.0, 1.0},
                                                                    {3, 0.0, 0.0, 0.0},
                                                                    {4, 0.5, 0.5, 0.0},
                                                                    {5, 0.0, 0.5, 0.5},
                                                                    {6, 0.0, 0.0, 0.5},
                                                                    {7, 0.5, 0.0, 0.0},
                                                                    {8, 0.5, 0.0, 0.5},
                                                                    {9, 0.0, 0.5, 0.0}}};

    matrix N_matrix(numerical_quadrature->points(), nodes());
    matrix local_quadrature_coordinates = matrix::Ones(numerical_quadrature->points(), 4);

    numerical_quadrature->evaluate([&](auto const& coordinate) {
        auto const& [l, r, s, t] = coordinate;

        auto const u = 1.0 - r - s - t;

        vector N(10);
        matrix rhea(10, 3);

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

        local_quadrature_coordinates(l, 0) = r;
        local_quadrature_coordinates(l, 1) = s;
        local_quadrature_coordinates(l, 2) = t;

        N_matrix.row(l) = N;

        return std::make_tuple(N, rhea);
    });

    // Compute extrapolation algorithm matrices
    matrix local_nodal_coordinates = matrix::Ones(nodes(), 4);

    for (auto const& [a, r, s, t] : local_coordinates)
    {
        local_nodal_coordinates(a, 0) = r;
        local_nodal_coordinates(a, 1) = s;
        local_nodal_coordinates(a, 2) = t;
    }
    compute_extrapolation_matrix(N_matrix, local_nodal_coordinates, local_quadrature_coordinates);
}
}
