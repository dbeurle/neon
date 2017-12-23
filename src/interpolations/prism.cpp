
#include "prism.hpp"

#include <array>
#include <memory>
#include <tuple>

namespace neon
{
prism6::prism6(PrismQuadrature::Rule rule)
    : VolumeInterpolation(std::make_unique<PrismQuadrature>(rule))
{
    this->precompute_shape_functions();
}

void prism6::precompute_shape_functions()
{
    using coordinate_type = std::tuple<int, double, double, double>;

    // Initialize nodal coordinates array as Xi, Eta, Zeta
    std::array<coordinate_type, 6> constexpr local_coordinates{{{0, 1.0, 0.0, -1.0},
                                                                {1, 0.0, 1.0, -1.0},
                                                                {2, 0.0, 0.0, -1.0},
                                                                {3, 1.0, 0.0, 1.0},
                                                                {4, 0.0, 1.0, 1.0},
                                                                {5, 0.0, 0.0, 1.0}}};

    matrix N_matrix(numerical_quadrature->points(), nodes());
    matrix local_quadrature_coordinates = matrix::Ones(numerical_quadrature->points(), 4);

    numerical_quadrature->evaluate([&](auto const& coordinate) {
        auto const & [l, r, s, xi] = coordinate;

        vector N(6);
        matrix rhea(6, 3);

        auto const t = 1.0 - r - s;

        N(1) = 0.5 * (1.0 - xi) * r;
        N(0) = 0.5 * (1.0 - xi) * s;
        N(2) = 0.5 * (1.0 - xi) * t;
        N(3) = 0.5 * (1.0 + xi) * r;
        N(4) = 0.5 * (1.0 + xi) * s;
        N(5) = 0.5 * (1.0 + xi) * t;

        rhea(0, 0) = 0.5 * (1.0 - xi);
        rhea(1, 0) = 0.0;
        rhea(2, 0) = -0.5 * (1.0 - xi);
        rhea(3, 0) = 0.5 * (1.0 + xi);
        rhea(4, 0) = 0.0;
        rhea(5, 0) = -0.5 * (1.0 + xi);

        rhea(0, 1) = 0.0;
        rhea(1, 1) = 0.5 * (1.0 - xi);
        rhea(2, 1) = -0.5 * (1.0 - xi);
        rhea(3, 1) = 0.0;
        rhea(4, 1) = 0.5 * (1.0 + xi);
        rhea(5, 1) = -0.5 * (1.0 + xi);

        rhea(0, 2) = -0.5 * r;
        rhea(1, 2) = -0.5 * s;
        rhea(2, 2) = -0.5 * t;
        rhea(3, 2) = 0.5 * r;
        rhea(4, 2) = 0.5 * s;
        rhea(5, 2) = 0.5 * t;

        local_quadrature_coordinates(l, 0) = r;
        local_quadrature_coordinates(l, 1) = s;
        local_quadrature_coordinates(l, 2) = xi;

        N_matrix.row(l) = N;

        return std::make_tuple(N, rhea);
    });

    // Compute extrapolation algorithm matrices
    matrix local_nodal_coordinates = matrix::Ones(nodes(), 4);

    for (auto const & [a, r_a, s_a, xi_a] : local_coordinates)
    {
        local_nodal_coordinates(a, 0) = r_a;
        local_nodal_coordinates(a, 1) = s_a;
        local_nodal_coordinates(a, 2) = xi_a;
    }
    compute_extrapolation_matrix(N_matrix, local_nodal_coordinates, local_quadrature_coordinates);
}

double prism6::compute_measure(matrix3x const& nodal_coordinates) const
{
    return numerical_quadrature->integrate(0.0, [&](auto const& femval, auto const& l) {
        auto const & [N, dN] = femval;

        matrix3 const Jacobian = nodal_coordinates * dN;

        return Jacobian.determinant();
    });
}
}
