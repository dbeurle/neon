
#include "UnitSphereQuadrature.hpp"

namespace neon
{
UnitSphereQuadrature::UnitSphereQuadrature()
{
    // 21 point unit sphere Gaussian quadrature scheme, see Bazant and Oh (1986)
    // assuming orthogonal symmetry

    // Weightings for unit sphere
    constexpr double w1 = 0.0265214244093 * 2.0, w2 = 0.0199301476312 * 2.0,
                     w3 = 0.0250712367487 * 2.0;
    w = {w1, w1, w1, w2, w2, w2, w2, w2, w2, w3, w3, w3, w3, w3, w3, w3, w3, w3, w3, w3, w3};

    // Directional cosines from unit sphere integration
    constexpr double dc0 = 0.707106781187, dc1 = 0.387907304067, dc2 = 0.836095596749;
    clist = {{0, 1.0, 0.0, 0.0},    {1, 0.0, 1.0, 0.0},  {2, 0.0, 0.0, 1.0},   {3, dc0, dc0, 0.0},
             {4, dc0, -dc0, 0.0},   {5, dc0, 0.0, dc0},  {6, dc0, 0.0, -dc0},  {7, 0.0, dc0, dc0},
             {8, 0.0, dc0, -dc0},   {9, dc1, dc1, dc2},  {10, dc1, dc1, -dc2}, {11, dc1, -dc1, dc2},
             {12, dc1, -dc1, -dc2}, {13, dc1, dc2, dc1}, {14, dc1, dc2, -dc1}, {15, dc1, -dc2, dc1},
             {16, dc1, -dc2, -dc1}, {17, dc2, dc1, dc1}, {18, dc2, dc1, -dc1}, {19, dc2, -dc1, dc1},
             {20, dc2, -dc1, -dc1}};

    this->precompute_coordinates();
}

void UnitSphereQuadrature::precompute_coordinates()
{
    this->evaluate([](auto const& coordinate) {

        auto const & [ l, r1, r2, r3 ] = coordinate;

        Vector3 t(r1, r2, r3);

        return std::make_tuple(t, t * t.transpose());
    });
}
}
