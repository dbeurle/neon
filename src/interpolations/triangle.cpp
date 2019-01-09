
#include "triangle.hpp"

namespace neon
{
triangle3::triangle3() : surface_interpolation(3, 0, 0) {}

auto triangle3::evaluate(coordinate_type const& coordinate) const noexcept(false) -> value_type
{
    [[maybe_unused]] auto const& [l, r, s] = coordinate;

    vector N(3);
    matrix rhea(3, 2);

    N(0) = r;
    N(1) = s;
    N(2) = 1.0 - r - s;

    rhea(0, 0) = 1.0;
    rhea(1, 0) = 0.0;
    rhea(2, 0) = -1.0;

    rhea(0, 1) = 0.0;
    rhea(1, 1) = 1.0;
    rhea(2, 1) = -1.0;

    return {N, rhea};
}

triangle6::triangle6() : surface_interpolation(6, 0, 0)
{
    m_local_coordinates =
        {{0, 1.0, 0.0}, {1, 0.0, 1.0}, {2, 0.0, 0.0}, {3, 0.5, 0.5}, {4, 0.0, 0.5}, {5, 0.5, 0.0}};
}

auto triangle6::evaluate(coordinate_type const& coordinate) const noexcept(false) -> value_type
{
    [[maybe_unused]] auto const& [l, r, s] = coordinate;

    auto const t = 1.0 - r - s;

    vector N(6);
    matrix rhea(6, 2);

    N(0) = r * (2.0 * r - 1.0);
    N(1) = s * (2.0 * s - 1.0);
    N(2) = t * (2.0 * t - 1.0);
    N(3) = 4.0 * r * s;
    N(4) = 4.0 * s * t;
    N(5) = 4.0 * r * t;

    // r coordinates
    rhea(0, 0) = 4.0 * r - 1.0;
    rhea(1, 0) = 0.0;
    rhea(2, 0) = -4.0 * t + 1.0;
    rhea(3, 0) = 4.0 * s;
    rhea(4, 0) = -4.0 * s;
    rhea(5, 0) = 4.0 * t - 4.0 * r;

    // s coordinates
    rhea(0, 1) = 0.0;
    rhea(1, 1) = 4.0 * s - 1.0;
    rhea(2, 1) = -4.0 * t + 1.0;
    rhea(3, 1) = 4.0 * r;
    rhea(4, 1) = 4.0 * t - 4.0 * s;
    rhea(5, 1) = -4.0 * r;

    return {N, rhea};
}

}
