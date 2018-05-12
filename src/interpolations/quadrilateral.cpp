
#include "quadrilateral.hpp"

#include "geometry/projection.hpp"

#include <array>
#include <tuple>

namespace neon
{
quadrilateral4::quadrilateral4(quadrilateral_quadrature::Rule rule)
    : surface_interpolation(std::make_unique<quadrilateral_quadrature>(rule))
{
    this->precompute_shape_functions();
}

void quadrilateral4::precompute_shape_functions()
{
    using NodalCoordinate = std::tuple<int, double, double>;

    std::array<NodalCoordinate, 4> constexpr local_coordinates = {
        {{0, -1.0, -1.0}, {1, 1.0, -1.0}, {2, 1.0, 1.0}, {3, -1.0, 1.0}}};

    matrix N_matrix(numerical_quadrature->points(), nodes());
    matrix local_quadrature_coordinates = matrix::Ones(numerical_quadrature->points(), 3);

    numerical_quadrature->evaluate([&](auto const& coordinate) {
        auto const & [ l, xi, eta ] = coordinate;

        vector N(4);
        matrix rhea(4, 2);

        for (auto const & [ a, xi_a, eta_a ] : local_coordinates)
        {
            N(a) = 0.25 * (1.0 + xi_a * xi) * (1.0 + eta_a * eta);
            rhea(a, 0) = 0.25 * (1.0 + eta_a * eta) * xi_a;
            rhea(a, 1) = 0.25 * (1.0 + xi_a * xi) * eta_a;
        }

        local_quadrature_coordinates(l, 0) = xi;
        local_quadrature_coordinates(l, 1) = eta;

        N_matrix.row(l) = N;

        return std::make_tuple(N, rhea);
    });

    // Compute extrapolation algorithm matrices
    matrix local_nodal_coordinates = matrix::Ones(nodes(), 3);

    for (auto const & [ a, xi_a, eta_a ] : local_coordinates)
    {
        local_nodal_coordinates(a, 0) = xi_a;
        local_nodal_coordinates(a, 1) = eta_a;
    }
    compute_extrapolation_matrix(N_matrix, local_nodal_coordinates, local_quadrature_coordinates);
}

double quadrilateral4::compute_measure(matrix const& nodal_coordinates) const
{
    return numerical_quadrature->integrate(0.0, [&](auto const& femval, auto const& l) {
        auto const& [N, dN] = femval;

        matrix2 const Jacobian = geometry::project_to_plane(nodal_coordinates) * dN;

        return Jacobian.determinant();
    });
}

quadrilateral8::quadrilateral8(quadrilateral_quadrature::Rule rule)
    : surface_interpolation(std::make_unique<quadrilateral_quadrature>(rule))
{
    this->precompute_shape_functions();
}

void quadrilateral8::precompute_shape_functions()
{
    using NodalCoordinate = std::tuple<int, double, double>;

    std::array<NodalCoordinate, 8> constexpr local_coordinates = {{{0, -1.0, -1.0},
                                                                   {1, 1.0, -1.0},
                                                                   {2, 1.0, 1.0},
                                                                   {3, -1.0, 1.0},
                                                                   {4, 0.0, -1.0},
                                                                   {5, -1.0, 0.0},
                                                                   {6, 0.0, 1.0},
                                                                   {7, -1.0, 0.0}}};

    matrix N_matrix(numerical_quadrature->points(), nodes());
    matrix local_quadrature_coordinates = matrix::Ones(numerical_quadrature->points(), 3);

    numerical_quadrature->evaluate([&](auto const& coordinate) {
        auto const & [ l, xi, eta ] = coordinate;

        vector N(8);
        matrix rhea(8, 2);

        N(0) = 0.25 * (1.0 - xi) * (1.0 - eta) * (-1.0 - xi - eta);
        N(1) = 0.25 * (1.0 + xi) * (1.0 - eta) * (-1.0 + xi - eta);
        N(2) = 0.25 * (1.0 + xi) * (1.0 + eta) * (-1.0 + xi + eta);
        N(3) = 0.25 * (1.0 - xi) * (1.0 + eta) * (-1.0 - xi + eta);
        N(4) = 0.5 * (1.0 - xi * xi) * (1.0 - eta);
        N(5) = 0.5 * (1.0 - eta * eta) * (1.0 + xi);
        N(6) = 0.5 * (1.0 - xi * xi) * (1.0 + eta);
        N(7) = 0.5 * (1.0 - eta * eta) * (1.0 - xi);

        rhea(0, 0) = 0.25 * (1.0 - eta) * (2.0 * xi + eta);
        rhea(1, 0) = 0.25 * (1.0 - eta) * (2.0 * xi - eta);
        rhea(2, 0) = 0.25 * (1.0 + eta) * (2.0 * xi + eta);
        rhea(3, 0) = 0.25 * (1.0 + eta) * (2.0 * xi - eta);
        rhea(4, 0) = -xi * (1.0 - eta);
        rhea(5, 0) = 0.5 * (1.0 - eta * eta);
        rhea(6, 0) = -xi * (1.0 + eta);
        rhea(7, 0) = -0.5 * (1.0 - eta * eta);

        rhea(0, 1) = 0.25 * (1.0 - xi) * (2.0 * eta + xi);
        rhea(1, 1) = 0.25 * (1.0 + xi) * (2.0 * eta - xi);
        rhea(2, 1) = 0.25 * (1.0 + xi) * (2.0 * eta + xi);
        rhea(3, 1) = 0.25 * (1.0 - xi) * (2.0 * eta - xi);
        rhea(4, 1) = -0.5 * (1.0 - xi * xi);
        rhea(5, 1) = -eta * (1.0 + xi);
        rhea(6, 1) = 0.5 * (1.0 - xi * xi);
        rhea(7, 1) = -eta * (1.0 - xi);

        local_quadrature_coordinates(l, 0) = xi;
        local_quadrature_coordinates(l, 1) = eta;

        N_matrix.row(l) = N;

        return std::make_tuple(N, rhea);
    });

    // Compute extrapolation algorithm matrices
    matrix local_nodal_coordinates = matrix::Ones(nodes(), 3);

    for (auto const & [ a, xi_a, eta_a ] : local_coordinates)
    {
        local_nodal_coordinates(a, 0) = xi_a;
        local_nodal_coordinates(a, 1) = eta_a;
    }
    compute_extrapolation_matrix(N_matrix, local_nodal_coordinates, local_quadrature_coordinates);
}

double quadrilateral8::compute_measure(matrix const& nodal_coordinates) const
{
    return numerical_quadrature->integrate(0.0, [&](auto const& femval, auto const& l) {
        auto const & [ N, dN ] = femval;

        matrix2 const Jacobian = geometry::project_to_plane(nodal_coordinates) * dN;

        return Jacobian.determinant();
    });
}

quadrilateral9::quadrilateral9(quadrilateral_quadrature::Rule rule)
    : surface_interpolation(std::make_unique<quadrilateral_quadrature>(rule))
{
    this->precompute_shape_functions();
}

void quadrilateral9::precompute_shape_functions()
{
    using NodalCoordinate = std::tuple<int, double, double>;

    std::array<NodalCoordinate, 9> constexpr local_coordinates = {{{0, -1.0, -1.0},
                                                                   {1, 1.0, -1.0},
                                                                   {2, 1.0, 1.0},
                                                                   {3, -1.0, 1.0},
                                                                   {4, 0.0, -1.0},
                                                                   {5, -1.0, 0.0},
                                                                   {6, 0.0, 1.0},
                                                                   {7, -1.0, 0.0},
                                                                   {8, 0.0, 0.0}}};

    matrix N_matrix(numerical_quadrature->points(), nodes());
    matrix local_quadrature_coordinates = matrix::Ones(numerical_quadrature->points(), 3);

    numerical_quadrature->evaluate([&](auto const& coordinate) {
        auto const & [ l, xi, eta ] = coordinate;

        vector N(9);
        matrix rhea(9, 2);

        N(0) = eta * xi * (eta - 1) * (xi - 1) / 4.0;
        N(1) = eta * xi * (eta - 1) * (xi + 1) / 4.0;
        N(2) = eta * xi * (eta + 1) * (xi + 1) / 4.0;
        N(3) = eta * xi * (eta + 1) * (xi - 1) / 4.0;
        N(4) = -eta * (eta - 1) * (xi - 1) * (xi + 1) / 2.0;
        N(5) = -xi * (eta - 1) * (eta + 1) * (xi + 1) / 2.0;
        N(6) = -eta * (eta + 1) * (xi - 1) * (xi + 1) / 2.0;
        N(7) = -xi * (eta - 1) * (eta + 1) * (xi - 1) / 2.0;
        N(8) = (eta - 1) * (eta + 1) * (xi - 1) * (xi + 1);

        rhea(0, 0) = eta * (eta - 1) * (2 * xi - 1) / 4.0;
        rhea(1, 0) = eta * (eta - 1) * (2 * xi + 1) / 4.0;
        rhea(2, 0) = eta * (eta + 1) * (2 * xi + 1) / 4.0;
        rhea(3, 0) = eta * (eta + 1) * (2 * xi - 1) / 4.0;
        rhea(4, 0) = -eta * xi * (eta - 1);
        rhea(5, 0) = -(eta - 1) * (eta + 1) * (2 * xi + 1) / 2.0;
        rhea(6, 0) = -eta * xi * (eta + 1);
        rhea(7, 0) = -(eta - 1) * (eta + 1) * (2 * xi - 1) / 2.0;
        rhea(8, 0) = 2 * xi * (eta - 1) * (eta + 1);

        rhea(0, 1) = xi * (2 * eta - 1) * (xi - 1) / 4.0;
        rhea(1, 1) = xi * (2 * eta - 1) * (xi + 1) / 4.0;
        rhea(2, 1) = xi * (2 * eta + 1) * (xi + 1) / 4.0;
        rhea(3, 1) = xi * (2 * eta + 1) * (xi - 1) / 4.0;
        rhea(4, 1) = -(2 * eta - 1) * (xi - 1) * (xi + 1) / 2.0;
        rhea(5, 1) = -eta * xi * (xi + 1);
        rhea(6, 1) = -(2 * eta + 1) * (xi - 1) * (xi + 1) / 2.0;
        rhea(7, 1) = -eta * xi * (xi - 1);
        rhea(8, 1) = 2 * eta * (xi - 1) * (xi + 1);

        local_quadrature_coordinates(l, 0) = xi;
        local_quadrature_coordinates(l, 1) = eta;

        N_matrix.row(l) = N;

        return std::make_tuple(N, rhea);
    });

    // Compute extrapolation algorithm matrices
    matrix local_nodal_coordinates = matrix::Ones(nodes(), 3);

    for (auto const & [ a, xi_a, eta_a ] : local_coordinates)
    {
        local_nodal_coordinates(a, 0) = xi_a;
        local_nodal_coordinates(a, 1) = eta_a;
    }
    compute_extrapolation_matrix(N_matrix, local_nodal_coordinates, local_quadrature_coordinates);
}

double quadrilateral9::compute_measure(matrix const& nodal_coordinates) const
{
    return numerical_quadrature->integrate(0.0, [&](auto const& femval, auto const& l) {
        auto const & [ N, dN ] = femval;

        matrix2 const Jacobian = geometry::project_to_plane(nodal_coordinates) * dN;

        return Jacobian.determinant();
    });
}
}
