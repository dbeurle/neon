
#include "tetrahedron_quadrature.hpp"

#include <algorithm>

namespace neon
{
tetrahedron_quadrature::tetrahedron_quadrature(point const p)
{
    switch (p)
    {
        case point::one:
        {
            w = {1.0};
            clist = {{0, 0.25, 0.25, 0.25}};
            break;
        }
        case point::four:
        {
            w = {0.25, 0.25, 0.25, 0.25};
            constexpr auto c0 = 0.585410196624969;
            constexpr auto c1 = 0.138196601125011;
            clist = {{0, c0, c1, c1}, {1, c1, c0, c1}, {2, c1, c1, c0}, {3, c1, c1, c1}};
            break;
        }
        case point::five:
        {
            w = {-4.0 / 5.0, 9.0 / 20.0, 9.0 / 20.0, 9.0 / 20.0, 9.0 / 20.0};
            clist = {{0, 0.25, 0.25, 0.25},
                     {1, 0.5, 1.0 / 6.0, 1.0 / 6.0},
                     {2, 1.0 / 6.0, 0.5, 1.0 / 6.0},
                     {3, 1.0 / 6.0, 1.0 / 6.0, 0.5},
                     {4, 1.0 / 6.0, 1.0 / 6.0, 1.0 / 6.0}};
            break;
        }
    }

    // Convert the weightings to proper quadrature format
    std::transform(begin(w), end(w), begin(w), [](auto const i) { return i / 6.0; });
}
}
