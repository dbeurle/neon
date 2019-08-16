
#include "line.hpp"

#include <Eigen/Geometry>

namespace neon
{
line2::line2() : line_interpolation(2, 1, 1) { m_local_coordinates = {{0, -1.0}, {1, 1.0}}; }

auto line2::evaluate(coordinate_type const& coordinate) const noexcept(false) -> value_type
{
    [[maybe_unused]] auto const [index, xi] = coordinate;

    vector N(2);
    matrix dN(2, 1);

    N(0) = 1.0 / 2.0 * (1.0 - xi);
    N(1) = 1.0 / 2.0 * (1.0 + xi);

    dN(0, 0) = -1.0 / 2.0;
    dN(1, 0) = 1.0 / 2.0;

    return {N, dN};
}

line3::line3() : line_interpolation(3, 2, 2)
{
    m_local_coordinates = {{0, -1.0}, {1, 0.0}, {2, 1.0}};
}

auto line3::evaluate(coordinate_type const& coordinate) const noexcept(false) -> value_type
{
    [[maybe_unused]] auto const [index, xi] = coordinate;

    vector N(3);
    matrix dN(3, 1);

    N(0) = 1.0 / 2.0 * xi * (xi - 1.0);
    N(1) = 1.0 - std::pow(xi, 2);
    N(2) = 1.0 / 2.0 * xi * (xi + 1.0);

    dN(0, 0) = 1.0 / 2.0 * (2.0 * xi - 1.0);
    dN(1, 0) = -2.0 * xi;
    dN(2, 0) = 1.0 / 2.0 * (2.0 * xi + 1.0);

    return {N, dN};
}
}
