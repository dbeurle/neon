
#include "UnitSphereQuadrature.hpp"

namespace neon
{
UnitSphereQuadrature::UnitSphereQuadrature(Rule const rule)
{
    switch (rule)
    {
        case Rule::BO21:
        {
            // 21 point unit sphere Gaussian quadrature scheme assuming orthogonal
            // symmetry see Bazant and Oh (1986) Table 1 on Page 43

            // Weightings for unit sphere
            constexpr double w1 = 0.0265214244093 * 2.0, w2 = 0.0199301476312 * 2.0,
                             w3 = 0.0250712367487 * 2.0;
            w = {w1, w1, w1, w2, w2, w2, w2, w2, w2, w3, w3,
                 w3, w3, w3, w3, w3, w3, w3, w3, w3, w3};

            // Directional cosines from unit sphere integration
            constexpr double dc0 = 0.707106781187, dc1 = 0.387907304067, dc2 = 0.836095596749;
            clist = {{0, 1.0, 0.0, 0.0},    {1, 0.0, 1.0, 0.0},    {2, 0.0, 0.0, 1.0},
                     {3, dc0, dc0, 0.0},    {4, dc0, -dc0, 0.0},   {5, dc0, 0.0, dc0},
                     {6, dc0, 0.0, -dc0},   {7, 0.0, dc0, dc0},    {8, 0.0, dc0, -dc0},
                     {9, dc1, dc1, dc2},    {10, dc1, dc1, -dc2},  {11, dc1, -dc1, dc2},
                     {12, dc1, -dc1, -dc2}, {13, dc1, dc2, dc1},   {14, dc1, dc2, -dc1},
                     {15, dc1, -dc2, dc1},  {16, dc1, -dc2, -dc1}, {17, dc2, dc1, dc1},
                     {18, dc2, dc1, -dc1},  {19, dc2, -dc1, dc1},  {20, dc2, -dc1, -dc1}};

            break;
        }
        case Rule::BO33:
        {
            // 21 point unit sphere Gaussian quadrature scheme assuming orthogonal
            // symmetry see Bazant and Oh (1986) Table 1 on Page 43

            // Weightings for unit sphere
            constexpr double w1 = 2.0 * 0.00985353993433, w2 = 2.0 * 0.0162969685886,
                             w3 = 2.0 * 0.0134788844008, w4 = 2.0 * 0.0175759129880;

            w = {w1, w1, w1,                                     //
                 w2, w2, w2, w2, w2, w2,                         //
                 w3, w3, w3, w3, w3, w3, w3, w3, w3, w3, w3, w3, //
                 w4, w4, w4, w4, w4, w4, w4, w4, w4, w4, w4, w4};

            // Directional cosines from unit sphere integration
            constexpr double dc0 = 0.707106781187, dc1 = 0.933898956394, dc2 = 0.357537045978,
                             dc3 = 0.437263676092, dc4 = 0.785875915868;
            clist = {{0, 1.0, 0.0, 0.0}, // 1
                     {1, 0.0, 1.0, 0.0}, // 2
                     {2, 0.0, 0.0, 1.0}, // 3

                     {3, dc0, dc0, 0.0},  // 4
                     {4, dc0, -dc0, 0.0}, // 5
                     {5, dc0, 0.0, dc0},  // 6
                     {6, dc0, 0.0, -dc0}, // 7
                     {7, 0.0, dc0, dc0},  // 8
                     {8, 0.0, dc0, -dc0}, // 9

                     {9, dc1, dc2, 0.0},   // 10
                     {10, dc1, -dc2, 0.0}, // 11
                     {11, dc2, dc1, 0.0},  // 12
                     {12, dc2, -dc1, 0.0}, // 13

                     {13, dc1, 0.0, dc2},  // 14
                     {14, dc1, 0.0, -dc2}, // 15
                     {15, dc2, 0.0, dc1},  // 16
                     {16, dc2, 0.0, -dc1}, // 17

                     {17, 0.0, dc1, dc2},  // 18
                     {18, 0.0, dc1, -dc2}, // 19
                     {19, 0.0, dc2, dc1},  // 20
                     {20, 0.0, dc2, -dc1}, // 21

                     {21, dc3, dc3, dc4},   // 22
                     {22, dc3, dc3, -dc4},  // 23
                     {23, dc3, -dc3, dc4},  // 24
                     {24, dc3, -dc3, -dc4}, // 25

                     {25, dc3, dc4, dc3},   // 26
                     {26, dc3, dc4, -dc3},  // 27
                     {27, dc3, -dc4, dc3},  // 28
                     {28, dc3, -dc4, -dc3}, // 29

                     {29, dc4, dc3, dc3},    // 30
                     {30, dc4, dc3, -dc3},   // 31
                     {31, dc4, -dc3, dc3},   // 32
                     {32, dc4, -dc3, -dc3}}; // 33

            break;
        }
    }
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
