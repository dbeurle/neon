
#include "Hexahedron.hpp"

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
        auto const& [l, xi, eta, zeta] = coordinate;

        Vector N(8);
        Matrix rhea(8, 3);

        for (auto const& [a, xi_a, eta_a, zeta_a] : local_coordinates)
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

    for (auto const& [a, xi_a, eta_a, zeta_a] : local_coordinates)
    {
        local_nodal_coordinates(a, 0) = xi_a;
        local_nodal_coordinates(a, 1) = eta_a;
        local_nodal_coordinates(a, 2) = zeta_a;
    }
    compute_extrapolation_matrix(N_matrix, local_nodal_coordinates, local_quadrature_coordinates);
}

double Hexahedron8::compute_measure(Matrix const& nodal_coordinates) const
{
    return numerical_quadrature->integrate(0.0, [&](auto const& femval, auto const& l) {
        auto const& [N, dN] = femval;

        Matrix3 const Jacobian = nodal_coordinates * dN;

        return Jacobian.determinant();
    });
}

Hexahedron20::Hexahedron20(HexahedronQuadrature::Rule rule)
    : VolumeInterpolation(std::make_unique<HexahedronQuadrature>(rule))
{
    this->precompute_shape_functions();
}

void Hexahedron20::precompute_shape_functions()
{
    using NodalCoordinate = std::tuple<int, double, double, double>;

    // Initialize nodal coordinates array as Xi, Eta, Zeta
    constexpr std::array<NodalCoordinate, 20> local_coordinates = {
        {{0, -1.0, -1.0, -1.0}, {1, 1.0, -1.0, -1.0}, {2, 1.0, 1.0, -1.0},  {3, -1.0, 1.0, -1.0},
         {4, -1.0, -1.0, 1.0},  {5, 1.0, -1.0, 1.0},  {6, 1.0, 1.0, 1.0},   {7, -1.0, 1.0, 1.0},
         {8, 0.0, -1.0, -1.0},  {9, 1.0, 0.0, -1.0},  {10, 0.0, 1.0, -1.0}, {11, -1.0, 0.0, -1.0},
         {12, 0.0, -1.0, 1.0},  {13, 1.0, 0.0, 1.0},  {14, 0.0, 1.0, 1.0},  {15, -1.0, 0.0, 1.0},
         {16, -1.0, -1.0, 0.0}, {17, 1.0, -1.0, 0.0}, {18, 1.0, 1.0, 0.0},  {19, -1.0, 1.0, 0.0}}};

    Matrix N_matrix(numerical_quadrature->points(), nodes());
    Matrix local_quadrature_coordinates = Matrix::Ones(numerical_quadrature->points(), 4);

    numerical_quadrature->evaluate([&](auto const& coordinate) {
        auto const& [l, xi, eta, zeta] = coordinate;

        Vector N(20);
        Matrix rhea(20, 3);

        N(0) = (eta - 1) * (xi - 1) * (zeta - 1) * (eta + xi + zeta + 2) / 8.0;
        N(1) = -(eta - 1) * (xi + 1) * (zeta - 1) * (eta - xi + zeta + 2) / 8.0;
        N(2) = -(eta + 1) * (xi + 1) * (zeta - 1) * (eta + xi - zeta - 2) / 8.0;
        N(3) = (eta + 1) * (xi - 1) * (zeta - 1) * (eta - xi - zeta - 2) / 8.0;
        N(4) = -(eta - 1) * (xi - 1) * (zeta + 1) * (eta + xi - zeta + 2) / 8.0;
        N(5) = (eta - 1) * (xi + 1) * (zeta + 1) * (eta - xi - zeta + 2) / 8.0;
        N(6) = (eta + 1) * (xi + 1) * (zeta + 1) * (eta + xi + zeta - 2) / 8.0;
        N(7) = -(eta + 1) * (xi - 1) * (zeta + 1) * (eta - xi + zeta - 2) / 8.0;
        N(8) = -(eta - 1) * (xi - 1) * (xi + 1) * (zeta - 1) / 4.0;
        N(9) = (eta - 1) * (eta + 1) * (xi + 1) * (zeta - 1) / 4.0;
        N(10) = (eta + 1) * (xi - 1) * (xi + 1) * (zeta - 1) / 4.0;
        N(11) = -(eta - 1) * (eta + 1) * (xi - 1) * (zeta - 1) / 4.0;
        N(12) = (eta - 1) * (xi - 1) * (xi + 1) * (zeta + 1) / 4.0;
        N(13) = -(eta - 1) * (eta + 1) * (xi + 1) * (zeta + 1) / 4.0;
        N(14) = -(eta + 1) * (xi - 1) * (xi + 1) * (zeta + 1) / 4.0;
        N(15) = (eta - 1) * (eta + 1) * (xi - 1) * (zeta + 1) / 4.0;
        N(16) = -(eta - 1) * (xi - 1) * (zeta - 1) * (zeta + 1) / 4.0;
        N(17) = (eta - 1) * (xi + 1) * (zeta - 1) * (zeta + 1) / 4.0;
        N(18) = -(eta + 1) * (xi + 1) * (zeta - 1) * (zeta + 1) / 4.0;
        N(19) = (eta + 1) * (xi - 1) * (zeta - 1) * (zeta + 1) / 4.0;

        rhea(0, 0) = (eta - 1) * (zeta - 1) * (eta + 2 * xi + zeta + 1) / 8.0;
        rhea(1, 0) = -(eta - 1) * (zeta - 1) * (eta - 2 * xi + zeta + 1) / 8.0;
        rhea(2, 0) = -(eta + 1) * (zeta - 1) * (eta + 2 * xi - zeta - 1) / 8.0;
        rhea(3, 0) = (eta + 1) * (zeta - 1) * (eta - 2 * xi - zeta - 1) / 8.0;
        rhea(4, 0) = -(eta - 1) * (zeta + 1) * (eta + 2 * xi - zeta + 1) / 8.0;
        rhea(5, 0) = (eta - 1) * (zeta + 1) * (eta - 2 * xi - zeta + 1) / 8.0;
        rhea(6, 0) = (eta + 1) * (zeta + 1) * (eta + 2 * xi + zeta - 1) / 8.0;
        rhea(7, 0) = -(eta + 1) * (zeta + 1) * (eta - 2 * xi + zeta - 1) / 8.0;
        rhea(8, 0) = -xi * (eta - 1) * (zeta - 1) / 2.0;
        rhea(9, 0) = (eta - 1) * (eta + 1) * (zeta - 1) / 4.0;
        rhea(10, 0) = xi * (eta + 1) * (zeta - 1) / 2.0;
        rhea(11, 0) = -(eta - 1) * (eta + 1) * (zeta - 1) / 4.0;
        rhea(12, 0) = xi * (eta - 1) * (zeta + 1) / 2.0;
        rhea(13, 0) = -(eta - 1) * (eta + 1) * (zeta + 1) / 4.0;
        rhea(14, 0) = -xi * (eta + 1) * (zeta + 1) / 2.0;
        rhea(15, 0) = (eta - 1) * (eta + 1) * (zeta + 1) / 4.0;
        rhea(16, 0) = -(eta - 1) * (zeta - 1) * (zeta + 1) / 4.0;
        rhea(17, 0) = (eta - 1) * (zeta - 1) * (zeta + 1) / 4.0;
        rhea(18, 0) = -(eta + 1) * (zeta - 1) * (zeta + 1) / 4.0;
        rhea(19, 0) = (eta + 1) * (zeta - 1) * (zeta + 1) / 4.0;

        rhea(0, 1) = (xi - 1) * (zeta - 1) * (2 * eta + xi + zeta + 1) / 8.0;
        rhea(1, 1) = -(xi + 1) * (zeta - 1) * (2 * eta - xi + zeta + 1) / 8.0;
        rhea(2, 1) = -(xi + 1) * (zeta - 1) * (2 * eta + xi - zeta - 1) / 8.0;
        rhea(3, 1) = (xi - 1) * (zeta - 1) * (2 * eta - xi - zeta - 1) / 8.0;
        rhea(4, 1) = -(xi - 1) * (zeta + 1) * (2 * eta + xi - zeta + 1) / 8.0;
        rhea(5, 1) = (xi + 1) * (zeta + 1) * (2 * eta - xi - zeta + 1) / 8.0;
        rhea(6, 1) = (xi + 1) * (zeta + 1) * (2 * eta + xi + zeta - 1) / 8.0;
        rhea(7, 1) = -(xi - 1) * (zeta + 1) * (2 * eta - xi + zeta - 1) / 8.0;
        rhea(8, 1) = -(xi - 1) * (xi + 1) * (zeta - 1) / 4.0;
        rhea(9, 1) = eta * (xi + 1) * (zeta - 1) / 2.0;
        rhea(10, 1) = (xi - 1) * (xi + 1) * (zeta - 1) / 4.0;
        rhea(11, 1) = -eta * (xi - 1) * (zeta - 1) / 2.0;
        rhea(12, 1) = (xi - 1) * (xi + 1) * (zeta + 1) / 4.0;
        rhea(13, 1) = -eta * (xi + 1) * (zeta + 1) / 2.0;
        rhea(14, 1) = -(xi - 1) * (xi + 1) * (zeta + 1) / 4.0;
        rhea(15, 1) = eta * (xi - 1) * (zeta + 1) / 2.0;
        rhea(16, 1) = -(xi - 1) * (zeta - 1) * (zeta + 1) / 4.0;
        rhea(17, 1) = (xi + 1) * (zeta - 1) * (zeta + 1) / 4.0;
        rhea(18, 1) = -(xi + 1) * (zeta - 1) * (zeta + 1) / 4.0;
        rhea(19, 1) = (xi - 1) * (zeta - 1) * (zeta + 1) / 4.0;

        rhea(0, 2) = (eta - 1) * (xi - 1) * (eta + xi + 2 * zeta + 1) / 8.0;
        rhea(1, 2) = -(eta - 1) * (xi + 1) * (eta - xi + 2 * zeta + 1) / 8.0;
        rhea(2, 2) = -(eta + 1) * (xi + 1) * (eta + xi - 2 * zeta - 1) / 8.0;
        rhea(3, 2) = (eta + 1) * (xi - 1) * (eta - xi - 2 * zeta - 1) / 8.0;
        rhea(4, 2) = -(eta - 1) * (xi - 1) * (eta + xi - 2 * zeta + 1) / 8.0;
        rhea(5, 2) = (eta - 1) * (xi + 1) * (eta - xi - 2 * zeta + 1) / 8.0;
        rhea(6, 2) = (eta + 1) * (xi + 1) * (eta + xi + 2 * zeta - 1) / 8.0;
        rhea(7, 2) = -(eta + 1) * (xi - 1) * (eta - xi + 2 * zeta - 1) / 8.0;
        rhea(8, 2) = -(eta - 1) * (xi - 1) * (xi + 1) / 4.0;
        rhea(9, 2) = (eta - 1) * (eta + 1) * (xi + 1) / 4.0;
        rhea(10, 2) = (eta + 1) * (xi - 1) * (xi + 1) / 4.0;
        rhea(11, 2) = -(eta - 1) * (eta + 1) * (xi - 1) / 4.0;
        rhea(12, 2) = (eta - 1) * (xi - 1) * (xi + 1) / 4.0;
        rhea(13, 2) = -(eta - 1) * (eta + 1) * (xi + 1) / 4.0;
        rhea(14, 2) = -(eta + 1) * (xi - 1) * (xi + 1) / 4.0;
        rhea(15, 2) = (eta - 1) * (eta + 1) * (xi - 1) / 4.0;
        rhea(16, 2) = -zeta * (eta - 1) * (xi - 1) / 2.0;
        rhea(17, 2) = zeta * (eta - 1) * (xi + 1) / 2.0;
        rhea(18, 2) = -zeta * (eta + 1) * (xi + 1) / 2.0;
        rhea(19, 2) = zeta * (eta + 1) * (xi - 1) / 2.0;

        local_quadrature_coordinates(l, 0) = xi;
        local_quadrature_coordinates(l, 1) = eta;
        local_quadrature_coordinates(l, 2) = zeta;

        N_matrix.row(l) = N;

        return std::make_tuple(N, rhea);
    });

    // Compute extrapolation algorithm matrices
    Matrix local_nodal_coordinates = Matrix::Ones(nodes(), 4);

    for (auto const& [a, xi_a, eta_a, zeta_a] : local_coordinates)
    {
        local_nodal_coordinates(a, 0) = xi_a;
        local_nodal_coordinates(a, 1) = eta_a;
        local_nodal_coordinates(a, 2) = zeta_a;
    }
    compute_extrapolation_matrix(N_matrix, local_nodal_coordinates, local_quadrature_coordinates);
}

double Hexahedron20::compute_measure(Matrix const& nodal_coordinates) const
{
    return numerical_quadrature->integrate(0.0, [&](auto const& femval, auto const& l) {
        auto const& [N, dN] = femval;

        Matrix3 const Jacobian = nodal_coordinates * dN;

        return Jacobian.determinant();
    });
}

Hexahedron27::Hexahedron27(HexahedronQuadrature::Rule rule)
    : VolumeInterpolation(std::make_unique<HexahedronQuadrature>(rule))
{
    this->precompute_shape_functions();
}

void Hexahedron27::precompute_shape_functions()
{
    using NodalCoordinate = std::tuple<int, double, double, double>;

    // Initialize nodal coordinates array as Xi, Eta, Zeta
    constexpr std::array<NodalCoordinate, 27> local_coordinates = {
        {{0, -1.0, -1.0, -1.0}, {1, 1.0, -1.0, -1.0}, {2, 1.0, 1.0, -1.0},  {3, -1.0, 1.0, -1.0},
         {4, -1.0, -1.0, 1.0},  {5, 1.0, -1.0, 1.0},  {6, 1.0, 1.0, 1.0},   {7, -1.0, 1.0, 1.0},
         {8, 0.0, -1.0, -1.0},  {9, 1.0, 0.0, -1.0},  {10, 0.0, 1.0, -1.0}, {11, -1.0, 0.0, -1.0},
         {12, 0.0, -1.0, 1.0},  {13, 1.0, 0.0, 1.0},  {14, 0.0, 1.0, 1.0},  {15, -1.0, 0.0, 1.0},
         {16, -1.0, -1.0, 0.0}, {17, 1.0, -1.0, 0.0}, {18, 1.0, 1.0, 0.0},  {19, -1.0, 1.0, 0.0},
         {20, 0.0, 0.0, -1.0},  {21, 0.0, 0.0, 1.0},  {22, 0.0, -1.0, 0.0}, {23, 0.0, 1.0, 0.0},
         {24, -1.0, 0.0, 0.0},  {25, 1.0, 0.0, 0.0},  {26, 0.0, 0.0, 0.0}}};

    Matrix N_matrix(numerical_quadrature->points(), nodes());
    Matrix local_quadrature_coordinates = Matrix::Ones(numerical_quadrature->points(), 4);

    numerical_quadrature->evaluate([&](auto const& coordinate) {
        auto const& [l, xi, eta, zeta] = coordinate;

        Vector N(27);
        Matrix rhea(27, 3);

        // Lagrange polynomial shape functions
        N(0) = eta * xi * zeta * (eta - 1) * (xi - 1) * (zeta - 1) / 8.0;
        N(1) = eta * xi * zeta * (eta - 1) * (xi + 1) * (zeta - 1) / 8.0;
        N(2) = eta * xi * zeta * (eta + 1) * (xi + 1) * (zeta - 1) / 8.0;
        N(3) = eta * xi * zeta * (eta + 1) * (xi - 1) * (zeta - 1) / 8.0;
        N(4) = eta * xi * zeta * (eta - 1) * (xi - 1) * (zeta + 1) / 8.0;
        N(5) = eta * xi * zeta * (eta - 1) * (xi + 1) * (zeta + 1) / 8.0;
        N(6) = eta * xi * zeta * (eta + 1) * (xi + 1) * (zeta + 1) / 8.0;
        N(7) = eta * xi * zeta * (eta + 1) * (xi - 1) * (zeta + 1) / 8.0;
        N(8) = -eta * zeta * (eta - 1) * (xi - 1) * (xi + 1) * (zeta - 1) / 4.0;
        N(9) = -xi * zeta * (eta - 1) * (eta + 1) * (xi + 1) * (zeta - 1) / 4.0;
        N(10) = -eta * zeta * (eta + 1) * (xi - 1) * (xi + 1) * (zeta - 1) / 4.0;
        N(11) = -xi * zeta * (eta - 1) * (eta + 1) * (xi - 1) * (zeta - 1) / 4.0;
        N(12) = -eta * zeta * (eta - 1) * (xi - 1) * (xi + 1) * (zeta + 1) / 4.0;
        N(13) = -xi * zeta * (eta - 1) * (eta + 1) * (xi + 1) * (zeta + 1) / 4.0;
        N(14) = -eta * zeta * (eta + 1) * (xi - 1) * (xi + 1) * (zeta + 1) / 4.0;
        N(15) = -xi * zeta * (eta - 1) * (eta + 1) * (xi - 1) * (zeta + 1) / 4.0;
        N(16) = -eta * xi * (eta - 1) * (xi - 1) * (zeta - 1) * (zeta + 1) / 4.0;
        N(17) = -eta * xi * (eta - 1) * (xi + 1) * (zeta - 1) * (zeta + 1) / 4.0;
        N(18) = -eta * xi * (eta + 1) * (xi + 1) * (zeta - 1) * (zeta + 1) / 4.0;
        N(19) = -eta * xi * (eta + 1) * (xi - 1) * (zeta - 1) * (zeta + 1) / 4.0;
        N(20) = zeta * (eta - 1) * (eta + 1) * (xi - 1) * (xi + 1) * (zeta - 1) / 2.0;
        N(21) = zeta * (eta - 1) * (eta + 1) * (xi - 1) * (xi + 1) * (zeta + 1) / 2.0;
        N(22) = eta * (eta - 1) * (xi - 1) * (xi + 1) * (zeta - 1) * (zeta + 1) / 2.0;
        N(23) = eta * (eta + 1) * (xi - 1) * (xi + 1) * (zeta - 1) * (zeta + 1) / 2.0;
        N(24) = xi * (eta - 1) * (eta + 1) * (xi - 1) * (zeta - 1) * (zeta + 1) / 2.0;
        N(25) = xi * (eta - 1) * (eta + 1) * (xi + 1) * (zeta - 1) * (zeta + 1) / 2.0;
        N(26) = -(eta - 1) * (eta + 1) * (xi - 1) * (xi + 1) * (zeta - 1) * (zeta + 1);

        // Lagrange polynomial shape function derivatives
        rhea(0, 0) = eta * zeta * (eta - 1) * (2 * xi - 1) * (zeta - 1) / 8.0;
        rhea(1, 0) = eta * zeta * (eta - 1) * (2 * xi + 1) * (zeta - 1) / 8.0;
        rhea(2, 0) = eta * zeta * (eta + 1) * (2 * xi + 1) * (zeta - 1) / 8.0;
        rhea(3, 0) = eta * zeta * (eta + 1) * (2 * xi - 1) * (zeta - 1) / 8.0;
        rhea(4, 0) = eta * zeta * (eta - 1) * (2 * xi - 1) * (zeta + 1) / 8.0;
        rhea(5, 0) = eta * zeta * (eta - 1) * (2 * xi + 1) * (zeta + 1) / 8.0;
        rhea(6, 0) = eta * zeta * (eta + 1) * (2 * xi + 1) * (zeta + 1) / 8.0;
        rhea(7, 0) = eta * zeta * (eta + 1) * (2 * xi - 1) * (zeta + 1) / 8.0;
        rhea(8, 0) = -eta * xi * zeta * (eta - 1) * (zeta - 1) / 2.0;
        rhea(9, 0) = -zeta * (eta - 1) * (eta + 1) * (2 * xi + 1) * (zeta - 1) / 4.0;
        rhea(10, 0) = -eta * xi * zeta * (eta + 1) * (zeta - 1) / 2.0;
        rhea(11, 0) = zeta * (eta - 1) * (eta + 1) * (-2 * xi + 1) * (zeta - 1) / 4.0;
        rhea(12, 0) = -eta * xi * zeta * (eta - 1) * (zeta + 1) / 2.0;
        rhea(13, 0) = -zeta * (eta - 1) * (eta + 1) * (2 * xi + 1) * (zeta + 1) / 4.0;
        rhea(14, 0) = -eta * xi * zeta * (eta + 1) * (zeta + 1) / 2.0;
        rhea(15, 0) = zeta * (eta - 1) * (eta + 1) * (-2 * xi + 1) * (zeta + 1) / 4.0;
        rhea(16, 0) = eta * (eta - 1) * (-2 * xi + 1) * (zeta - 1) * (zeta + 1) / 4.0;
        rhea(17, 0) = -eta * (eta - 1) * (2 * xi + 1) * (zeta - 1) * (zeta + 1) / 4.0;
        rhea(18, 0) = -eta * (eta + 1) * (2 * xi + 1) * (zeta - 1) * (zeta + 1) / 4.0;
        rhea(19, 0) = eta * (eta + 1) * (-2 * xi + 1) * (zeta - 1) * (zeta + 1) / 4.0;
        rhea(20, 0) = xi * zeta * (eta - 1) * (eta + 1) * (zeta - 1);
        rhea(21, 0) = xi * zeta * (eta - 1) * (eta + 1) * (zeta + 1);
        rhea(22, 0) = eta * xi * (eta - 1) * (zeta - 1) * (zeta + 1);
        rhea(23, 0) = eta * xi * (eta + 1) * (zeta - 1) * (zeta + 1);
        rhea(24, 0) = (eta - 1) * (eta + 1) * (2 * xi - 1) * (zeta - 1) * (zeta + 1) / 2.0;
        rhea(25, 0) = (eta - 1) * (eta + 1) * (2 * xi + 1) * (zeta - 1) * (zeta + 1) / 2.0;
        rhea(26, 0) = -2 * xi * (eta - 1) * (eta + 1) * (zeta - 1) * (zeta + 1);

        rhea(0, 1) = xi * zeta * (2 * eta - 1) * (xi - 1) * (zeta - 1) / 8.0;
        rhea(1, 1) = xi * zeta * (2 * eta - 1) * (xi + 1) * (zeta - 1) / 8.0;
        rhea(2, 1) = xi * zeta * (2 * eta + 1) * (xi + 1) * (zeta - 1) / 8.0;
        rhea(3, 1) = xi * zeta * (2 * eta + 1) * (xi - 1) * (zeta - 1) / 8.0;
        rhea(4, 1) = xi * zeta * (2 * eta - 1) * (xi - 1) * (zeta + 1) / 8.0;
        rhea(5, 1) = xi * zeta * (2 * eta - 1) * (xi + 1) * (zeta + 1) / 8.0;
        rhea(6, 1) = xi * zeta * (2 * eta + 1) * (xi + 1) * (zeta + 1) / 8.0;
        rhea(7, 1) = xi * zeta * (2 * eta + 1) * (xi - 1) * (zeta + 1) / 8.0;
        rhea(8, 1) = zeta * (-2 * eta + 1) * (xi - 1) * (xi + 1) * (zeta - 1) / 4.0;
        rhea(9, 1) = -eta * xi * zeta * (xi + 1) * (zeta - 1) / 2.0;
        rhea(10, 1) = -zeta * (2 * eta + 1) * (xi - 1) * (xi + 1) * (zeta - 1) / 4.0;
        rhea(11, 1) = -eta * xi * zeta * (xi - 1) * (zeta - 1) / 2.0;
        rhea(12, 1) = zeta * (-2 * eta + 1) * (xi - 1) * (xi + 1) * (zeta + 1) / 4.0;
        rhea(13, 1) = -eta * xi * zeta * (xi + 1) * (zeta + 1) / 2.0;
        rhea(14, 1) = -zeta * (2 * eta + 1) * (xi - 1) * (xi + 1) * (zeta + 1) / 4.0;
        rhea(15, 1) = -eta * xi * zeta * (xi - 1) * (zeta + 1) / 2.0;
        rhea(16, 1) = xi * (-2 * eta + 1) * (xi - 1) * (zeta - 1) * (zeta + 1) / 4.0;
        rhea(17, 1) = xi * (-2 * eta + 1) * (xi + 1) * (zeta - 1) * (zeta + 1) / 4.0;
        rhea(18, 1) = -xi * (2 * eta + 1) * (xi + 1) * (zeta - 1) * (zeta + 1) / 4.0;
        rhea(19, 1) = -xi * (2 * eta + 1) * (xi - 1) * (zeta - 1) * (zeta + 1) / 4.0;
        rhea(20, 1) = eta * zeta * (xi - 1) * (xi + 1) * (zeta - 1);
        rhea(21, 1) = eta * zeta * (xi - 1) * (xi + 1) * (zeta + 1);
        rhea(22, 1) = (2 * eta - 1) * (xi - 1) * (xi + 1) * (zeta - 1) * (zeta + 1) / 2.0;
        rhea(23, 1) = (2 * eta + 1) * (xi - 1) * (xi + 1) * (zeta - 1) * (zeta + 1) / 2.0;
        rhea(24, 1) = eta * xi * (xi - 1) * (zeta - 1) * (zeta + 1);
        rhea(25, 1) = eta * xi * (xi + 1) * (zeta - 1) * (zeta + 1);
        rhea(26, 1) = -2 * eta * (xi - 1) * (xi + 1) * (zeta - 1) * (zeta + 1);

        rhea(0, 2) = eta * xi * (eta - 1) * (xi - 1) * (2 * zeta - 1) / 8.0;
        rhea(1, 2) = eta * xi * (eta - 1) * (xi + 1) * (2 * zeta - 1) / 8.0;
        rhea(2, 2) = eta * xi * (eta + 1) * (xi + 1) * (2 * zeta - 1) / 8.0;
        rhea(3, 2) = eta * xi * (eta + 1) * (xi - 1) * (2 * zeta - 1) / 8.0;
        rhea(4, 2) = eta * xi * (eta - 1) * (xi - 1) * (2 * zeta + 1) / 8.0;
        rhea(5, 2) = eta * xi * (eta - 1) * (xi + 1) * (2 * zeta + 1) / 8.0;
        rhea(6, 2) = eta * xi * (eta + 1) * (xi + 1) * (2 * zeta + 1) / 8.0;
        rhea(7, 2) = eta * xi * (eta + 1) * (xi - 1) * (2 * zeta + 1) / 8.0;
        rhea(8, 2) = eta * (eta - 1) * (xi - 1) * (xi + 1) * (-2 * zeta + 1) / 4.0;
        rhea(9, 2) = xi * (eta - 1) * (eta + 1) * (xi + 1) * (-2 * zeta + 1) / 4.0;
        rhea(10, 2) = eta * (eta + 1) * (xi - 1) * (xi + 1) * (-2 * zeta + 1) / 4.0;
        rhea(11, 2) = xi * (eta - 1) * (eta + 1) * (xi - 1) * (-2 * zeta + 1) / 4.0;
        rhea(12, 2) = -eta * (eta - 1) * (xi - 1) * (xi + 1) * (2 * zeta + 1) / 4.0;
        rhea(13, 2) = -xi * (eta - 1) * (eta + 1) * (xi + 1) * (2 * zeta + 1) / 4.0;
        rhea(14, 2) = -eta * (eta + 1) * (xi - 1) * (xi + 1) * (2 * zeta + 1) / 4.0;
        rhea(15, 2) = -xi * (eta - 1) * (eta + 1) * (xi - 1) * (2 * zeta + 1) / 4.0;
        rhea(16, 2) = -eta * xi * zeta * (eta - 1) * (xi - 1) / 2.0;
        rhea(17, 2) = -eta * xi * zeta * (eta - 1) * (xi + 1) / 2.0;
        rhea(18, 2) = -eta * xi * zeta * (eta + 1) * (xi + 1) / 2.0;
        rhea(19, 2) = -eta * xi * zeta * (eta + 1) * (xi - 1) / 2.0;
        rhea(20, 2) = (eta - 1) * (eta + 1) * (xi - 1) * (xi + 1) * (2 * zeta - 1) / 2.0;
        rhea(21, 2) = (eta - 1) * (eta + 1) * (xi - 1) * (xi + 1) * (2 * zeta + 1) / 2.0;
        rhea(22, 2) = eta * zeta * (eta - 1) * (xi - 1) * (xi + 1);
        rhea(23, 2) = eta * zeta * (eta + 1) * (xi - 1) * (xi + 1);
        rhea(24, 2) = xi * zeta * (eta - 1) * (eta + 1) * (xi - 1);
        rhea(25, 2) = xi * zeta * (eta - 1) * (eta + 1) * (xi + 1);
        rhea(26, 2) = -2 * zeta * (eta - 1) * (eta + 1) * (xi - 1) * (xi + 1);

        local_quadrature_coordinates(l, 0) = xi;
        local_quadrature_coordinates(l, 1) = eta;
        local_quadrature_coordinates(l, 2) = zeta;

        N_matrix.row(l) = N;

        return std::make_tuple(N, rhea);
    });

    // Compute extrapolation algorithm matrices
    Matrix local_nodal_coordinates = Matrix::Ones(nodes(), 4);

    for (auto const& [a, xi_a, eta_a, zeta_a] : local_coordinates)
    {
        local_nodal_coordinates(a, 0) = xi_a;
        local_nodal_coordinates(a, 1) = eta_a;
        local_nodal_coordinates(a, 2) = zeta_a;
    }
    compute_extrapolation_matrix(N_matrix, local_nodal_coordinates, local_quadrature_coordinates);
}

double Hexahedron27::compute_measure(Matrix const& nodal_coordinates) const
{
    return numerical_quadrature->integrate(0.0, [&](auto const& femval, auto const& l) {
        auto const& [N, dN] = femval;

        Matrix3 const Jacobian = nodal_coordinates * dN;

        return Jacobian.determinant();
    });
}
}
