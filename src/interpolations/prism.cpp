
#include "prism.hpp"

#include <array>
#include <memory>
#include <tuple>

namespace neon
{
prism6::prism6() : volume_interpolation(6, 1, 1)
{
    m_local_coordinates = {{0, 1.0, 0.0, -1.0},
                           {1, 0.0, 1.0, -1.0},
                           {2, 0.0, 0.0, -1.0},
                           {3, 1.0, 0.0, 1.0},
                           {4, 0.0, 1.0, 1.0},
                           {5, 0.0, 0.0, 1.0}};
}

auto prism6::evaluate(coordinate_type const& coordinate) const noexcept(false) -> value_type
{
    [[maybe_unused]] auto const& [l, r, s, xi] = coordinate;

    vector N(6);
    matrix dN(6, 3);

    auto const t = 1.0 - r - s;

    N(1) = 0.5 * (1.0 - xi) * r;
    N(0) = 0.5 * (1.0 - xi) * s;
    N(2) = 0.5 * (1.0 - xi) * t;
    N(3) = 0.5 * (1.0 + xi) * r;
    N(4) = 0.5 * (1.0 + xi) * s;
    N(5) = 0.5 * (1.0 + xi) * t;

    dN(0, 0) = 0.5 * (1.0 - xi);
    dN(1, 0) = 0.0;
    dN(2, 0) = -0.5 * (1.0 - xi);
    dN(3, 0) = 0.5 * (1.0 + xi);
    dN(4, 0) = 0.0;
    dN(5, 0) = -0.5 * (1.0 + xi);

    dN(0, 1) = 0.0;
    dN(1, 1) = 0.5 * (1.0 - xi);
    dN(2, 1) = -0.5 * (1.0 - xi);
    dN(3, 1) = 0.0;
    dN(4, 1) = 0.5 * (1.0 + xi);
    dN(5, 1) = -0.5 * (1.0 + xi);

    dN(0, 2) = -0.5 * r;
    dN(1, 2) = -0.5 * s;
    dN(2, 2) = -0.5 * t;
    dN(3, 2) = 0.5 * r;
    dN(4, 2) = 0.5 * s;
    dN(5, 2) = 0.5 * t;

    return {N, dN};
}

prism15::prism15() : volume_interpolation(15, 3, 4)
{
    m_local_coordinates = {{0, 1.0, 0.0, -1.0},
                           {1, 0.0, 1.0, -1.0},
                           {2, 0.0, 0.0, -1.0},
                           {3, 0.5, 0.5, -1.0},
                           {4, 0.0, 0.5, -1.0},
                           {5, 0.5, 0.0, -1.0},
                           //
                           {6, 1.0, 0.0, 0.0},
                           {7, 0.0, 1.0, 0.0},
                           {8, 0.0, 0.0, 0.0},
                           //
                           {9, 1.0, 0.0, 1.0},
                           {10, 0.0, 1.0, 1.0},
                           {11, 0.0, 0.0, 1.0},
                           {12, 0.5, 0.5, 1.0},
                           {13, 0.0, 0.5, 1.0},
                           {14, 0.5, 0.0, 1.0}};
}

auto prism15::evaluate(coordinate_type const& coordinate) const noexcept(false) -> value_type
{
    [[maybe_unused]] auto const& [l, r, s, xi] = coordinate;

    vector N(15);
    matrix dN(15, 3);

    N(0) = r * xi * (2 * r - 1) * (xi - 1) / 2.0;
    N(1) = s * xi * (2 * s - 1) * (xi - 1) / 2.0;
    N(2) = xi * (xi - 1) * (r + s - 1) * (2 * r + 2 * s - 1) / 2.0;
    N(3) = 2 * r * s * xi * (xi - 1);
    N(4) = -2 * s * xi * (xi - 1) * (r + s - 1);
    N(5) = -2 * r * xi * (xi - 1) * (r + s - 1);
    N(6) = -r * (xi - 1) * (xi + 1);
    N(7) = -s * (xi - 1) * (xi + 1);
    N(8) = (xi - 1) * (xi + 1) * (r + s - 1);
    N(9) = r * xi * (2 * r - 1) * (xi + 1) / 2.0;
    N(10) = s * xi * (2 * s - 1) * (xi + 1) / 2.0;
    N(11) = xi * (xi + 1) * (r + s - 1) * (2 * r + 2 * s - 1) / 2.0;
    N(12) = 2 * r * s * xi * (xi + 1);
    N(13) = -2 * s * xi * (xi + 1) * (r + s - 1);
    N(14) = -2 * r * xi * (xi + 1) * (r + s - 1);

    dN(0, 0) = xi * (4 * r - 1) * (xi - 1) / 2.0;
    dN(1, 0) = 0;
    dN(2, 0) = xi * (xi - 1) * (4 * r + 4 * s - 3) / 2.0;
    dN(3, 0) = 2 * s * xi * (xi - 1);
    dN(4, 0) = -2 * s * xi * (xi - 1);
    dN(5, 0) = -2 * xi * (xi - 1) * (2 * r + s - 1);
    dN(6, 0) = -(xi - 1) * (xi + 1);
    dN(7, 0) = 0;
    dN(8, 0) = (xi - 1) * (xi + 1);
    dN(9, 0) = xi * (4 * r - 1) * (xi + 1) / 2.0;
    dN(10, 0) = 0;
    dN(11, 0) = xi * (xi + 1) * (4 * r + 4 * s - 3) / 2.0;
    dN(12, 0) = 2 * s * xi * (xi + 1);
    dN(13, 0) = -2 * s * xi * (xi + 1);
    dN(14, 0) = -2 * xi * (xi + 1) * (2 * r + s - 1);

    dN(0, 1) = 0;
    dN(1, 1) = xi * (4 * s - 1) * (xi - 1) / 2.0;
    dN(2, 1) = xi * (xi - 1) * (4 * r + 4 * s - 3) / 2.0;
    dN(3, 1) = 2 * r * xi * (xi - 1);
    dN(4, 1) = -2 * xi * (xi - 1) * (r + 2 * s - 1);
    dN(5, 1) = -2 * r * xi * (xi - 1);
    dN(6, 1) = 0;
    dN(7, 1) = -(xi - 1) * (xi + 1);
    dN(8, 1) = (xi - 1) * (xi + 1);
    dN(9, 1) = 0;
    dN(10, 1) = xi * (4 * s - 1) * (xi + 1) / 2.0;
    dN(11, 1) = xi * (xi + 1) * (4 * r + 4 * s - 3) / 2.0;
    dN(12, 1) = 2 * r * xi * (xi + 1);
    dN(13, 1) = -2 * xi * (xi + 1) * (r + 2 * s - 1);
    dN(14, 1) = -2 * r * xi * (xi + 1);

    dN(0, 2) = r * (2 * r - 1) * (2 * xi - 1) / 2.0;
    dN(1, 2) = s * (2 * s - 1) * (2 * xi - 1) / 2.0;
    dN(2, 2) = (2 * xi - 1) * (r + s - 1) * (2 * r + 2 * s - 1) / 2.0;
    dN(3, 2) = 2 * r * s * (2 * xi - 1);
    dN(4, 2) = -2 * s * (2 * xi - 1) * (r + s - 1);
    dN(5, 2) = -2 * r * (2 * xi - 1) * (r + s - 1);
    dN(6, 2) = -2 * r * xi;
    dN(7, 2) = -2 * s * xi;
    dN(8, 2) = 2 * xi * (r + s - 1);
    dN(9, 2) = r * (2 * r - 1) * (2 * xi + 1) / 2.0;
    dN(10, 2) = s * (2 * s - 1) * (2 * xi + 1) / 2.0;
    dN(11, 2) = (2 * xi + 1) * (r + s - 1) * (2 * r + 2 * s - 1) / 2.0;
    dN(12, 2) = 2 * r * s * (2 * xi + 1);
    dN(13, 2) = -2 * s * (2 * xi + 1) * (r + s - 1);
    dN(14, 2) = -2 * r * (2 * xi + 1) * (r + s - 1);

    return {N, dN};
}
}
