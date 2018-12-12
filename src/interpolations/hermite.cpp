
#include "hermite.hpp"

namespace neon
{
hermite::hermite(line_quadrature::point const p)
    : line_interpolation(std::make_unique<line_quadrature>(p))
{
    this->precompute_shape_functions();
}

void hermite::precompute_shape_functions()
{
    using coordinates_type = std::tuple<int, double>;

    // Initialize nodal coordinates array as Xi, Eta, Zeta
    std::array<coordinates_type, 2> constexpr local_coordinates{{{0, -1.0}, {1, 1.0}}};

    matrix N_matrix(m_quadrature->points(), nodes());
    matrix local_quadrature_coordinates = matrix::Ones(m_quadrature->points(), 2);

    m_quadrature->evaluate([&](auto const& coordinates) {
        auto const& [l, xi] = coordinates;

        vector N(4);
        matrix dN(4, 1);

        // shape function
        N(0) = 1.0 / 4.0 * std::pow(1.0 - xi, 2) * (xi + 2);
        N(1) = 1.0 / 4.0 * std::pow(1.0 + xi, 2) * (2 - xi);
        N(2) = 1.0 / 8.0 * std::pow(1.0 - xi, 2) * (xi + 1.0);
        N(3) = 1.0 / 8.0 * std::pow(1.0 + xi, 2) * (xi - 1.0);

        // FIXME second derivative of the shape functions
        dN(0, 0) = -1.0 / 2.0;
        dN(1, 0) = 1.0 / 2.0;
        dN(2, 0) = -1.0 / 2.0;
        dN(3, 0) = 1.0 / 2.0;

        local_quadrature_coordinates(l, 0) = xi;

        N_matrix.row(l) = N;

        return std::make_tuple(N, dN);
    });

    // Compute extrapolation algorithm matrices
    matrix local_nodal_coordinates = matrix::Ones(nodes(), 2);

    for (auto const& [a, xi_a] : local_coordinates)
    {
        local_nodal_coordinates(a, 0) = xi_a;
    }
    compute_extrapolation_matrix(N_matrix, local_nodal_coordinates, local_quadrature_coordinates);
}

double hermite::compute_measure(matrix const& nodal_coordinates) const
{
    return (nodal_coordinates.col(0) - nodal_coordinates.col(1)).norm();
}
}
