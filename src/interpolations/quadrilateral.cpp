
#include "quadrilateral.hpp"

#include "geometry/projection.hpp"

#include <array>
#include <tuple>

namespace neon
{
quadrilateral4::quadrilateral4() : surface_interpolation(4, 1, 1)
{
    m_local_coordinates = {{0, -1.0, -1.0}, {1, 1.0, -1.0}, {2, 1.0, 1.0}, {3, -1.0, 1.0}};
}

auto quadrilateral4::evaluate(coordinate_type const& coordinate) const noexcept(false) -> value_type
{
    auto const& [l, xi, eta] = coordinate;

    vector N(4);
    matrix dN(4, 2);

    for (auto const& [a, xi_a, eta_a] : m_local_coordinates)
    {
        N(a) = 0.25 * (1.0 + xi_a * xi) * (1.0 + eta_a * eta);
        dN(a, 0) = 0.25 * (1.0 + eta_a * eta) * xi_a;
        dN(a, 1) = 0.25 * (1.0 + xi_a * xi) * eta_a;
    }

    return {N, dN};
}

quadrilateral8::quadrilateral8() : surface_interpolation(8, 0, 0)
{
    m_local_coordinates = {{0, -1.0, -1.0},
                           {1, 1.0, -1.0},
                           {2, 1.0, 1.0},
                           {3, -1.0, 1.0},
                           {4, 0.0, -1.0},
                           {5, -1.0, 0.0},
                           {6, 0.0, 1.0},
                           {7, -1.0, 0.0}};
}

auto quadrilateral8::evaluate(coordinate_type const& coordinate) const noexcept(false) -> value_type
{
    [[maybe_unused]] auto const& [l, xi, eta] = coordinate;

    vector N(8);
    matrix dN(8, 2);

    N(0) = 0.25 * (1.0 - xi) * (1.0 - eta) * (-1.0 - xi - eta);
    N(1) = 0.25 * (1.0 + xi) * (1.0 - eta) * (-1.0 + xi - eta);
    N(2) = 0.25 * (1.0 + xi) * (1.0 + eta) * (-1.0 + xi + eta);
    N(3) = 0.25 * (1.0 - xi) * (1.0 + eta) * (-1.0 - xi + eta);
    N(4) = 0.5 * (1.0 - xi * xi) * (1.0 - eta);
    N(5) = 0.5 * (1.0 - eta * eta) * (1.0 + xi);
    N(6) = 0.5 * (1.0 - xi * xi) * (1.0 + eta);
    N(7) = 0.5 * (1.0 - eta * eta) * (1.0 - xi);

    dN(0, 0) = 0.25 * (1.0 - eta) * (2.0 * xi + eta);
    dN(1, 0) = 0.25 * (1.0 - eta) * (2.0 * xi - eta);
    dN(2, 0) = 0.25 * (1.0 + eta) * (2.0 * xi + eta);
    dN(3, 0) = 0.25 * (1.0 + eta) * (2.0 * xi - eta);
    dN(4, 0) = -xi * (1.0 - eta);
    dN(5, 0) = 0.5 * (1.0 - eta * eta);
    dN(6, 0) = -xi * (1.0 + eta);
    dN(7, 0) = -0.5 * (1.0 - eta * eta);

    dN(0, 1) = 0.25 * (1.0 - xi) * (2.0 * eta + xi);
    dN(1, 1) = 0.25 * (1.0 + xi) * (2.0 * eta - xi);
    dN(2, 1) = 0.25 * (1.0 + xi) * (2.0 * eta + xi);
    dN(3, 1) = 0.25 * (1.0 - xi) * (2.0 * eta - xi);
    dN(4, 1) = -0.5 * (1.0 - xi * xi);
    dN(5, 1) = -eta * (1.0 + xi);
    dN(6, 1) = 0.5 * (1.0 - xi * xi);
    dN(7, 1) = -eta * (1.0 - xi);

    return {N, dN};
}

quadrilateral9::quadrilateral9() : surface_interpolation(9, 0, 0)
{
    m_local_coordinates = {{0, -1.0, -1.0},
                           {1, 1.0, -1.0},
                           {2, 1.0, 1.0},
                           {3, -1.0, 1.0},
                           {4, 0.0, -1.0},
                           {5, -1.0, 0.0},
                           {6, 0.0, 1.0},
                           {7, -1.0, 0.0},
                           {8, 0.0, 0.0}};
}

auto quadrilateral9::evaluate(coordinate_type const& coordinate) const noexcept(false) -> value_type
{
    [[maybe_unused]] auto const& [l, xi, eta] = coordinate;

    vector N(9);
    matrix dN(9, 2);

    N(0) = eta * xi * (eta - 1) * (xi - 1) / 4.0;
    N(1) = eta * xi * (eta - 1) * (xi + 1) / 4.0;
    N(2) = eta * xi * (eta + 1) * (xi + 1) / 4.0;
    N(3) = eta * xi * (eta + 1) * (xi - 1) / 4.0;
    N(4) = -eta * (eta - 1) * (xi - 1) * (xi + 1) / 2.0;
    N(5) = -xi * (eta - 1) * (eta + 1) * (xi + 1) / 2.0;
    N(6) = -eta * (eta + 1) * (xi - 1) * (xi + 1) / 2.0;
    N(7) = -xi * (eta - 1) * (eta + 1) * (xi - 1) / 2.0;
    N(8) = (eta - 1) * (eta + 1) * (xi - 1) * (xi + 1);

    dN(0, 0) = eta * (eta - 1) * (2 * xi - 1) / 4.0;
    dN(1, 0) = eta * (eta - 1) * (2 * xi + 1) / 4.0;
    dN(2, 0) = eta * (eta + 1) * (2 * xi + 1) / 4.0;
    dN(3, 0) = eta * (eta + 1) * (2 * xi - 1) / 4.0;
    dN(4, 0) = -eta * xi * (eta - 1);
    dN(5, 0) = -(eta - 1) * (eta + 1) * (2 * xi + 1) / 2.0;
    dN(6, 0) = -eta * xi * (eta + 1);
    dN(7, 0) = -(eta - 1) * (eta + 1) * (2 * xi - 1) / 2.0;
    dN(8, 0) = 2 * xi * (eta - 1) * (eta + 1);

    dN(0, 1) = xi * (2 * eta - 1) * (xi - 1) / 4.0;
    dN(1, 1) = xi * (2 * eta - 1) * (xi + 1) / 4.0;
    dN(2, 1) = xi * (2 * eta + 1) * (xi + 1) / 4.0;
    dN(3, 1) = xi * (2 * eta + 1) * (xi - 1) / 4.0;
    dN(4, 1) = -(2 * eta - 1) * (xi - 1) * (xi + 1) / 2.0;
    dN(5, 1) = -eta * xi * (xi + 1);
    dN(6, 1) = -(2 * eta + 1) * (xi - 1) * (xi + 1) / 2.0;
    dN(7, 1) = -eta * xi * (xi - 1);
    dN(8, 1) = 2 * eta * (xi - 1) * (xi + 1);

    return {N, dN};
}
}
