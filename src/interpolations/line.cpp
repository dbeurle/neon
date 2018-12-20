
#include "line.hpp"

#include <Eigen/Geometry>

namespace neon
{
line2::line2(line_quadrature::point const p)
    : line_interpolation(std::make_unique<line_quadrature>(p), 2)
{
    this->precompute_shape_functions();
}

void line2::precompute_shape_functions()
{
    using coordinates_type = std::tuple<int, double>;

    // Initialize nodal coordinates array as Xi, Eta, Zeta
    std::array<coordinates_type, 2> constexpr local_coordinates{{{0, -1.0}, {1, 1.0}}};

    matrix N_matrix(m_quadrature->points(), number_of_nodes());
    matrix local_quadrature_coordinates = matrix::Ones(m_quadrature->points(), 2);

    m_quadrature->evaluate([&](auto const& coordinates) {
        auto const& [l, xi] = coordinates;

        vector N(2);
        matrix dN(2, 1);

        N(0) = 1.0 / 2.0 * (1.0 - xi);
        N(1) = 1.0 / 2.0 * (1.0 + xi);

        dN(0, 0) = -1.0 / 2.0;
        dN(1, 0) = 1.0 / 2.0;

        local_quadrature_coordinates(l, 0) = xi;

        N_matrix.row(l) = N;

        return std::make_tuple(N, dN);
    });

    // Compute extrapolation algorithm matrices
    matrix local_nodal_coordinates = matrix::Ones(number_of_nodes(), 2);

    for (auto const& [a, xi_a] : local_coordinates)
    {
        local_nodal_coordinates(a, 0) = xi_a;
    }
    compute_extrapolation_matrix(N_matrix, local_nodal_coordinates, local_quadrature_coordinates);
}

double line2::compute_measure(matrix const& nodal_coordinates) const
{
    return (nodal_coordinates.col(0) - nodal_coordinates.col(1)).norm();
}

line3::line3(line_quadrature::point const p)
    : line_interpolation(std::make_unique<line_quadrature>(p), 3)
{
    this->precompute_shape_functions();
}

void line3::precompute_shape_functions()
{
    using coordinates_type = std::tuple<int, double>;

    // Initialize nodal coordinates array as Xi, Eta, Zeta
    std::array<coordinates_type, 3> constexpr local_coordinates{{{0, -1.0}, {1, 0.0}, {2, 1.0}}};

    matrix N_matrix(m_quadrature->points(), number_of_nodes());
    matrix local_quadrature_coordinates = matrix::Ones(m_quadrature->points(), 2);

    m_quadrature->evaluate([&](auto const& coordinates) {
        auto const& [l, xi] = coordinates;

        vector N(3);
        matrix rhea(3, 1);

        N(0) = 1.0 / 2.0 * xi * (xi - 1.0);
        N(1) = 1.0 - std::pow(xi, 2);
        N(2) = 1.0 / 2.0 * xi * (xi + 1.0);

        rhea(0, 0) = 1.0 / 2.0 * (2.0 * xi - 1.0);
        rhea(1, 0) = -2.0 * xi;
        rhea(2, 0) = 1.0 / 2.0 * (2.0 * xi + 1.0);

        local_quadrature_coordinates(l, 0) = xi;

        N_matrix.row(l) = N;

        return std::make_tuple(N, rhea);
    });

    // Compute extrapolation algorithm matrices
    matrix local_nodal_coordinates = matrix::Ones(number_of_nodes(), 2);

    for (auto const& [a, xi_a] : local_coordinates)
    {
        local_nodal_coordinates(a, 0) = xi_a;
    }
    compute_extrapolation_matrix(N_matrix, local_nodal_coordinates, local_quadrature_coordinates);
}

double line3::compute_measure(matrix const& nodal_coordinates) const
{
    return (nodal_coordinates.col(0) - nodal_coordinates.col(2)).norm();
}
}
