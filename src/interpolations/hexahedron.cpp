
#include "hexahedron.hpp"

#include <array>
#include <memory>
#include <tuple>

namespace neon
{
hexahedron8::hexahedron8() : volume_interpolation(8, 1, 3)
{
    m_local_coordinates = {{0, -1.0, -1.0, -1.0},
                           {1, 1.0, -1.0, -1.0},
                           {2, 1.0, 1.0, -1.0},
                           {3, -1.0, 1.0, -1.0},
                           {4, -1.0, -1.0, 1.0},
                           {5, 1.0, -1.0, 1.0},
                           {6, 1.0, 1.0, 1.0},
                           {7, -1.0, 1.0, 1.0}};
}

auto hexahedron8::evaluate(coordinate_type const& coordinate) const noexcept(false) -> value_type
{
    [[maybe_unused]] auto const [index, xi, eta, zeta] = coordinate;

    vector N(8);
    matrix dN(8, 3);

    for (auto const& [a, xi_a, eta_a, zeta_a] : m_local_coordinates)
    {
        N(a) = 1.0 / 8.0 * (1.0 + xi_a * xi) * (1.0 + eta_a * eta) * (1.0 + zeta_a * zeta);

        dN(a, 0) = 1.0 / 8.0 * xi_a * (1.0 + eta_a * eta) * (1.0 + zeta_a * zeta);
        dN(a, 1) = 1.0 / 8.0 * (1.0 + xi_a * xi) * eta_a * (1.0 + zeta_a * zeta);
        dN(a, 2) = 1.0 / 8.0 * (1.0 + xi_a * xi) * (1.0 + eta_a * eta) * zeta_a;
    }
    return {N, dN};
}

hexahedron20::hexahedron20() : volume_interpolation(20, 2, 4)
{
    m_local_coordinates = {{0, -1.0, -1.0, -1.0}, {1, 1.0, -1.0, -1.0},  {2, 1.0, 1.0, -1.0},
                           {3, -1.0, 1.0, -1.0},  {4, -1.0, -1.0, 1.0},  {5, 1.0, -1.0, 1.0},
                           {6, 1.0, 1.0, 1.0},    {7, -1.0, 1.0, 1.0},   {8, 0.0, -1.0, -1.0},
                           {9, 1.0, 0.0, -1.0},   {10, 0.0, 1.0, -1.0},  {11, -1.0, 0.0, -1.0},
                           {12, 0.0, -1.0, 1.0},  {13, 1.0, 0.0, 1.0},   {14, 0.0, 1.0, 1.0},
                           {15, -1.0, 0.0, 1.0},  {16, -1.0, -1.0, 0.0}, {17, 1.0, -1.0, 0.0},
                           {18, 1.0, 1.0, 0.0},   {19, -1.0, 1.0, 0.0}};
}

auto hexahedron20::evaluate(coordinate_type const& coordinate) const noexcept(false) -> value_type
{
    [[maybe_unused]] auto const [l, xi, eta, zeta] = coordinate;

    vector N(20);
    matrix dN(20, 3);

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

    dN(0, 0) = (eta - 1) * (zeta - 1) * (eta + 2 * xi + zeta + 1) / 8.0;
    dN(1, 0) = -(eta - 1) * (zeta - 1) * (eta - 2 * xi + zeta + 1) / 8.0;
    dN(2, 0) = -(eta + 1) * (zeta - 1) * (eta + 2 * xi - zeta - 1) / 8.0;
    dN(3, 0) = (eta + 1) * (zeta - 1) * (eta - 2 * xi - zeta - 1) / 8.0;
    dN(4, 0) = -(eta - 1) * (zeta + 1) * (eta + 2 * xi - zeta + 1) / 8.0;
    dN(5, 0) = (eta - 1) * (zeta + 1) * (eta - 2 * xi - zeta + 1) / 8.0;
    dN(6, 0) = (eta + 1) * (zeta + 1) * (eta + 2 * xi + zeta - 1) / 8.0;
    dN(7, 0) = -(eta + 1) * (zeta + 1) * (eta - 2 * xi + zeta - 1) / 8.0;
    dN(8, 0) = -xi * (eta - 1) * (zeta - 1) / 2.0;
    dN(9, 0) = (eta - 1) * (eta + 1) * (zeta - 1) / 4.0;
    dN(10, 0) = xi * (eta + 1) * (zeta - 1) / 2.0;
    dN(11, 0) = -(eta - 1) * (eta + 1) * (zeta - 1) / 4.0;
    dN(12, 0) = xi * (eta - 1) * (zeta + 1) / 2.0;
    dN(13, 0) = -(eta - 1) * (eta + 1) * (zeta + 1) / 4.0;
    dN(14, 0) = -xi * (eta + 1) * (zeta + 1) / 2.0;
    dN(15, 0) = (eta - 1) * (eta + 1) * (zeta + 1) / 4.0;
    dN(16, 0) = -(eta - 1) * (zeta - 1) * (zeta + 1) / 4.0;
    dN(17, 0) = (eta - 1) * (zeta - 1) * (zeta + 1) / 4.0;
    dN(18, 0) = -(eta + 1) * (zeta - 1) * (zeta + 1) / 4.0;
    dN(19, 0) = (eta + 1) * (zeta - 1) * (zeta + 1) / 4.0;

    dN(0, 1) = (xi - 1) * (zeta - 1) * (2 * eta + xi + zeta + 1) / 8.0;
    dN(1, 1) = -(xi + 1) * (zeta - 1) * (2 * eta - xi + zeta + 1) / 8.0;
    dN(2, 1) = -(xi + 1) * (zeta - 1) * (2 * eta + xi - zeta - 1) / 8.0;
    dN(3, 1) = (xi - 1) * (zeta - 1) * (2 * eta - xi - zeta - 1) / 8.0;
    dN(4, 1) = -(xi - 1) * (zeta + 1) * (2 * eta + xi - zeta + 1) / 8.0;
    dN(5, 1) = (xi + 1) * (zeta + 1) * (2 * eta - xi - zeta + 1) / 8.0;
    dN(6, 1) = (xi + 1) * (zeta + 1) * (2 * eta + xi + zeta - 1) / 8.0;
    dN(7, 1) = -(xi - 1) * (zeta + 1) * (2 * eta - xi + zeta - 1) / 8.0;
    dN(8, 1) = -(xi - 1) * (xi + 1) * (zeta - 1) / 4.0;
    dN(9, 1) = eta * (xi + 1) * (zeta - 1) / 2.0;
    dN(10, 1) = (xi - 1) * (xi + 1) * (zeta - 1) / 4.0;
    dN(11, 1) = -eta * (xi - 1) * (zeta - 1) / 2.0;
    dN(12, 1) = (xi - 1) * (xi + 1) * (zeta + 1) / 4.0;
    dN(13, 1) = -eta * (xi + 1) * (zeta + 1) / 2.0;
    dN(14, 1) = -(xi - 1) * (xi + 1) * (zeta + 1) / 4.0;
    dN(15, 1) = eta * (xi - 1) * (zeta + 1) / 2.0;
    dN(16, 1) = -(xi - 1) * (zeta - 1) * (zeta + 1) / 4.0;
    dN(17, 1) = (xi + 1) * (zeta - 1) * (zeta + 1) / 4.0;
    dN(18, 1) = -(xi + 1) * (zeta - 1) * (zeta + 1) / 4.0;
    dN(19, 1) = (xi - 1) * (zeta - 1) * (zeta + 1) / 4.0;

    dN(0, 2) = (eta - 1) * (xi - 1) * (eta + xi + 2 * zeta + 1) / 8.0;
    dN(1, 2) = -(eta - 1) * (xi + 1) * (eta - xi + 2 * zeta + 1) / 8.0;
    dN(2, 2) = -(eta + 1) * (xi + 1) * (eta + xi - 2 * zeta - 1) / 8.0;
    dN(3, 2) = (eta + 1) * (xi - 1) * (eta - xi - 2 * zeta - 1) / 8.0;
    dN(4, 2) = -(eta - 1) * (xi - 1) * (eta + xi - 2 * zeta + 1) / 8.0;
    dN(5, 2) = (eta - 1) * (xi + 1) * (eta - xi - 2 * zeta + 1) / 8.0;
    dN(6, 2) = (eta + 1) * (xi + 1) * (eta + xi + 2 * zeta - 1) / 8.0;
    dN(7, 2) = -(eta + 1) * (xi - 1) * (eta - xi + 2 * zeta - 1) / 8.0;
    dN(8, 2) = -(eta - 1) * (xi - 1) * (xi + 1) / 4.0;
    dN(9, 2) = (eta - 1) * (eta + 1) * (xi + 1) / 4.0;
    dN(10, 2) = (eta + 1) * (xi - 1) * (xi + 1) / 4.0;
    dN(11, 2) = -(eta - 1) * (eta + 1) * (xi - 1) / 4.0;
    dN(12, 2) = (eta - 1) * (xi - 1) * (xi + 1) / 4.0;
    dN(13, 2) = -(eta - 1) * (eta + 1) * (xi + 1) / 4.0;
    dN(14, 2) = -(eta + 1) * (xi - 1) * (xi + 1) / 4.0;
    dN(15, 2) = (eta - 1) * (eta + 1) * (xi - 1) / 4.0;
    dN(16, 2) = -zeta * (eta - 1) * (xi - 1) / 2.0;
    dN(17, 2) = zeta * (eta - 1) * (xi + 1) / 2.0;
    dN(18, 2) = -zeta * (eta + 1) * (xi + 1) / 2.0;
    dN(19, 2) = zeta * (eta + 1) * (xi - 1) / 2.0;

    return {N, dN};
}

hexahedron27::hexahedron27() : volume_interpolation(27, 4, 4)
{
    m_local_coordinates = {{0, -1.0, -1.0, -1.0}, {1, 1.0, -1.0, -1.0},  {2, 1.0, 1.0, -1.0},
                           {3, -1.0, 1.0, -1.0},  {4, -1.0, -1.0, 1.0},  {5, 1.0, -1.0, 1.0},
                           {6, 1.0, 1.0, 1.0},    {7, -1.0, 1.0, 1.0},   {8, 0.0, -1.0, -1.0},
                           {9, 1.0, 0.0, -1.0},   {10, 0.0, 1.0, -1.0},  {11, -1.0, 0.0, -1.0},
                           {12, 0.0, -1.0, 1.0},  {13, 1.0, 0.0, 1.0},   {14, 0.0, 1.0, 1.0},
                           {15, -1.0, 0.0, 1.0},  {16, -1.0, -1.0, 0.0}, {17, 1.0, -1.0, 0.0},
                           {18, 1.0, 1.0, 0.0},   {19, -1.0, 1.0, 0.0},  {20, 0.0, 0.0, -1.0},
                           {21, 0.0, 0.0, 1.0},   {22, 0.0, -1.0, 0.0},  {23, 0.0, 1.0, 0.0},
                           {24, -1.0, 0.0, 0.0},  {25, 1.0, 0.0, 0.0},   {26, 0.0, 0.0, 0.0}};
}

auto hexahedron27::evaluate(coordinate_type const& coordinate) const noexcept(false) -> value_type
{
    [[maybe_unused]] auto const& [l, xi, eta, zeta] = coordinate;

    vector N(27);
    matrix dN(27, 3);

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
    dN(0, 0) = eta * zeta * (eta - 1) * (2 * xi - 1) * (zeta - 1) / 8.0;
    dN(1, 0) = eta * zeta * (eta - 1) * (2 * xi + 1) * (zeta - 1) / 8.0;
    dN(2, 0) = eta * zeta * (eta + 1) * (2 * xi + 1) * (zeta - 1) / 8.0;
    dN(3, 0) = eta * zeta * (eta + 1) * (2 * xi - 1) * (zeta - 1) / 8.0;
    dN(4, 0) = eta * zeta * (eta - 1) * (2 * xi - 1) * (zeta + 1) / 8.0;
    dN(5, 0) = eta * zeta * (eta - 1) * (2 * xi + 1) * (zeta + 1) / 8.0;
    dN(6, 0) = eta * zeta * (eta + 1) * (2 * xi + 1) * (zeta + 1) / 8.0;
    dN(7, 0) = eta * zeta * (eta + 1) * (2 * xi - 1) * (zeta + 1) / 8.0;
    dN(8, 0) = -eta * xi * zeta * (eta - 1) * (zeta - 1) / 2.0;
    dN(9, 0) = -zeta * (eta - 1) * (eta + 1) * (2 * xi + 1) * (zeta - 1) / 4.0;
    dN(10, 0) = -eta * xi * zeta * (eta + 1) * (zeta - 1) / 2.0;
    dN(11, 0) = zeta * (eta - 1) * (eta + 1) * (-2 * xi + 1) * (zeta - 1) / 4.0;
    dN(12, 0) = -eta * xi * zeta * (eta - 1) * (zeta + 1) / 2.0;
    dN(13, 0) = -zeta * (eta - 1) * (eta + 1) * (2 * xi + 1) * (zeta + 1) / 4.0;
    dN(14, 0) = -eta * xi * zeta * (eta + 1) * (zeta + 1) / 2.0;
    dN(15, 0) = zeta * (eta - 1) * (eta + 1) * (-2 * xi + 1) * (zeta + 1) / 4.0;
    dN(16, 0) = eta * (eta - 1) * (-2 * xi + 1) * (zeta - 1) * (zeta + 1) / 4.0;
    dN(17, 0) = -eta * (eta - 1) * (2 * xi + 1) * (zeta - 1) * (zeta + 1) / 4.0;
    dN(18, 0) = -eta * (eta + 1) * (2 * xi + 1) * (zeta - 1) * (zeta + 1) / 4.0;
    dN(19, 0) = eta * (eta + 1) * (-2 * xi + 1) * (zeta - 1) * (zeta + 1) / 4.0;
    dN(20, 0) = xi * zeta * (eta - 1) * (eta + 1) * (zeta - 1);
    dN(21, 0) = xi * zeta * (eta - 1) * (eta + 1) * (zeta + 1);
    dN(22, 0) = eta * xi * (eta - 1) * (zeta - 1) * (zeta + 1);
    dN(23, 0) = eta * xi * (eta + 1) * (zeta - 1) * (zeta + 1);
    dN(24, 0) = (eta - 1) * (eta + 1) * (2 * xi - 1) * (zeta - 1) * (zeta + 1) / 2.0;
    dN(25, 0) = (eta - 1) * (eta + 1) * (2 * xi + 1) * (zeta - 1) * (zeta + 1) / 2.0;
    dN(26, 0) = -2 * xi * (eta - 1) * (eta + 1) * (zeta - 1) * (zeta + 1);

    dN(0, 1) = xi * zeta * (2 * eta - 1) * (xi - 1) * (zeta - 1) / 8.0;
    dN(1, 1) = xi * zeta * (2 * eta - 1) * (xi + 1) * (zeta - 1) / 8.0;
    dN(2, 1) = xi * zeta * (2 * eta + 1) * (xi + 1) * (zeta - 1) / 8.0;
    dN(3, 1) = xi * zeta * (2 * eta + 1) * (xi - 1) * (zeta - 1) / 8.0;
    dN(4, 1) = xi * zeta * (2 * eta - 1) * (xi - 1) * (zeta + 1) / 8.0;
    dN(5, 1) = xi * zeta * (2 * eta - 1) * (xi + 1) * (zeta + 1) / 8.0;
    dN(6, 1) = xi * zeta * (2 * eta + 1) * (xi + 1) * (zeta + 1) / 8.0;
    dN(7, 1) = xi * zeta * (2 * eta + 1) * (xi - 1) * (zeta + 1) / 8.0;
    dN(8, 1) = zeta * (-2 * eta + 1) * (xi - 1) * (xi + 1) * (zeta - 1) / 4.0;
    dN(9, 1) = -eta * xi * zeta * (xi + 1) * (zeta - 1) / 2.0;
    dN(10, 1) = -zeta * (2 * eta + 1) * (xi - 1) * (xi + 1) * (zeta - 1) / 4.0;
    dN(11, 1) = -eta * xi * zeta * (xi - 1) * (zeta - 1) / 2.0;
    dN(12, 1) = zeta * (-2 * eta + 1) * (xi - 1) * (xi + 1) * (zeta + 1) / 4.0;
    dN(13, 1) = -eta * xi * zeta * (xi + 1) * (zeta + 1) / 2.0;
    dN(14, 1) = -zeta * (2 * eta + 1) * (xi - 1) * (xi + 1) * (zeta + 1) / 4.0;
    dN(15, 1) = -eta * xi * zeta * (xi - 1) * (zeta + 1) / 2.0;
    dN(16, 1) = xi * (-2 * eta + 1) * (xi - 1) * (zeta - 1) * (zeta + 1) / 4.0;
    dN(17, 1) = xi * (-2 * eta + 1) * (xi + 1) * (zeta - 1) * (zeta + 1) / 4.0;
    dN(18, 1) = -xi * (2 * eta + 1) * (xi + 1) * (zeta - 1) * (zeta + 1) / 4.0;
    dN(19, 1) = -xi * (2 * eta + 1) * (xi - 1) * (zeta - 1) * (zeta + 1) / 4.0;
    dN(20, 1) = eta * zeta * (xi - 1) * (xi + 1) * (zeta - 1);
    dN(21, 1) = eta * zeta * (xi - 1) * (xi + 1) * (zeta + 1);
    dN(22, 1) = (2 * eta - 1) * (xi - 1) * (xi + 1) * (zeta - 1) * (zeta + 1) / 2.0;
    dN(23, 1) = (2 * eta + 1) * (xi - 1) * (xi + 1) * (zeta - 1) * (zeta + 1) / 2.0;
    dN(24, 1) = eta * xi * (xi - 1) * (zeta - 1) * (zeta + 1);
    dN(25, 1) = eta * xi * (xi + 1) * (zeta - 1) * (zeta + 1);
    dN(26, 1) = -2 * eta * (xi - 1) * (xi + 1) * (zeta - 1) * (zeta + 1);

    dN(0, 2) = eta * xi * (eta - 1) * (xi - 1) * (2 * zeta - 1) / 8.0;
    dN(1, 2) = eta * xi * (eta - 1) * (xi + 1) * (2 * zeta - 1) / 8.0;
    dN(2, 2) = eta * xi * (eta + 1) * (xi + 1) * (2 * zeta - 1) / 8.0;
    dN(3, 2) = eta * xi * (eta + 1) * (xi - 1) * (2 * zeta - 1) / 8.0;
    dN(4, 2) = eta * xi * (eta - 1) * (xi - 1) * (2 * zeta + 1) / 8.0;
    dN(5, 2) = eta * xi * (eta - 1) * (xi + 1) * (2 * zeta + 1) / 8.0;
    dN(6, 2) = eta * xi * (eta + 1) * (xi + 1) * (2 * zeta + 1) / 8.0;
    dN(7, 2) = eta * xi * (eta + 1) * (xi - 1) * (2 * zeta + 1) / 8.0;
    dN(8, 2) = eta * (eta - 1) * (xi - 1) * (xi + 1) * (-2 * zeta + 1) / 4.0;
    dN(9, 2) = xi * (eta - 1) * (eta + 1) * (xi + 1) * (-2 * zeta + 1) / 4.0;
    dN(10, 2) = eta * (eta + 1) * (xi - 1) * (xi + 1) * (-2 * zeta + 1) / 4.0;
    dN(11, 2) = xi * (eta - 1) * (eta + 1) * (xi - 1) * (-2 * zeta + 1) / 4.0;
    dN(12, 2) = -eta * (eta - 1) * (xi - 1) * (xi + 1) * (2 * zeta + 1) / 4.0;
    dN(13, 2) = -xi * (eta - 1) * (eta + 1) * (xi + 1) * (2 * zeta + 1) / 4.0;
    dN(14, 2) = -eta * (eta + 1) * (xi - 1) * (xi + 1) * (2 * zeta + 1) / 4.0;
    dN(15, 2) = -xi * (eta - 1) * (eta + 1) * (xi - 1) * (2 * zeta + 1) / 4.0;
    dN(16, 2) = -eta * xi * zeta * (eta - 1) * (xi - 1) / 2.0;
    dN(17, 2) = -eta * xi * zeta * (eta - 1) * (xi + 1) / 2.0;
    dN(18, 2) = -eta * xi * zeta * (eta + 1) * (xi + 1) / 2.0;
    dN(19, 2) = -eta * xi * zeta * (eta + 1) * (xi - 1) / 2.0;
    dN(20, 2) = (eta - 1) * (eta + 1) * (xi - 1) * (xi + 1) * (2 * zeta - 1) / 2.0;
    dN(21, 2) = (eta - 1) * (eta + 1) * (xi - 1) * (xi + 1) * (2 * zeta + 1) / 2.0;
    dN(22, 2) = eta * zeta * (eta - 1) * (xi - 1) * (xi + 1);
    dN(23, 2) = eta * zeta * (eta + 1) * (xi - 1) * (xi + 1);
    dN(24, 2) = xi * zeta * (eta - 1) * (eta + 1) * (xi - 1);
    dN(25, 2) = xi * zeta * (eta - 1) * (eta + 1) * (xi + 1);
    dN(26, 2) = -2 * zeta * (eta - 1) * (eta + 1) * (xi - 1) * (xi + 1);

    return {N, dN};
}
}
