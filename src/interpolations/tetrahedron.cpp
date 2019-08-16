
#include "tetrahedron.hpp"

namespace neon
{
tetrahedron4::tetrahedron4() : volume_interpolation(4, 1, 1) {}

auto tetrahedron4::evaluate(coordinate_type const& coordinate) const noexcept(false) -> value_type
{
    [[maybe_unused]] auto const& [l, r, s, t] = coordinate;

    vector N(4);
    matrix dN(4, 3);

    N(0) = r;
    N(1) = s;
    N(2) = t;
    N(3) = 1.0 - r - s - t;

    dN(0, 0) = 1.0;
    dN(0, 1) = 0.0;
    dN(0, 2) = 0.0;

    dN(1, 0) = 0.0;
    dN(1, 1) = 1.0;
    dN(1, 2) = 0.0;

    dN(2, 0) = 0.0;
    dN(2, 1) = 0.0;
    dN(2, 2) = 1.0;

    dN(3, 0) = -1.0;
    dN(3, 1) = -1.0;
    dN(3, 2) = -1.0;

    return {N, dN};
}

tetrahedron10::tetrahedron10() : volume_interpolation(10, 2, 2)
{
    m_local_coordinates = {{0, 1.0, 0.0, 0.0},
                           {1, 0.0, 1.0, 0.0},
                           {2, 0.0, 0.0, 1.0},
                           {3, 0.0, 0.0, 0.0},
                           {4, 0.5, 0.5, 0.0},
                           {5, 0.0, 0.5, 0.5},
                           {6, 0.0, 0.0, 0.5},
                           {7, 0.5, 0.0, 0.0},
                           {8, 0.5, 0.0, 0.5},
                           {9, 0.0, 0.5, 0.0}};
}

auto tetrahedron10::evaluate(coordinate_type const& coordinate) const noexcept(false) -> value_type
{
    [[maybe_unused]] auto const& [l, r, s, t] = coordinate;

    auto const u = 1.0 - r - s - t;

    vector N(10);
    matrix dN(10, 3);

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

    dN(0, 0) = 4.0 * r - 1.0;
    dN(1, 0) = 0.0;
    dN(2, 0) = 0.0;
    dN(3, 0) = -3.0 + 4.0 * r + 4.0 * s + 4.0 * t;
    dN(4, 0) = 4.0 * s;
    dN(5, 0) = 0.0;
    dN(6, 0) = -4.0 * t;
    dN(7, 0) = 4.0 - 4.0 * s - 8.0 * r - 4.0 * t;
    dN(8, 0) = 4.0 * t;
    dN(9, 0) = -4.0 * s;

    dN(0, 1) = 0.0;
    dN(1, 1) = 4.0 * s - 1.0;
    dN(2, 1) = 0.0;
    dN(3, 1) = -3.0 + 4.0 * r + 4.0 * s + 4.0 * t;
    dN(4, 1) = 4.0 * r;
    dN(5, 1) = 4.0 * t;
    dN(6, 1) = -4.0 * t;
    dN(7, 1) = -4.0 * r;
    dN(8, 1) = 0.0;
    dN(9, 1) = 4.0 - 4.0 * r - 8.0 * s - 4.0 * t;

    dN(0, 2) = 0.0;
    dN(1, 2) = 0.0;
    dN(2, 2) = 4.0 * t - 1.0;
    dN(3, 2) = -3.0 + 4.0 * r + 4.0 * s + 4.0 * t;
    dN(4, 2) = 0.0;
    dN(5, 2) = 4.0 * s;
    dN(6, 2) = 4.0 - 4.0 * r - 4.0 * s - 8.0 * t;
    dN(7, 2) = -4.0 * r;
    dN(8, 2) = 4.0 * r;
    dN(9, 2) = -4.0 * s;

    return {N, dN};
}
}
