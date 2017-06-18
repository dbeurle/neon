
#include "TetrahedronQuadrature.hpp"

#include <range/v3/action.hpp>

namespace neon
{
TetrahedronQuadrature::TetrahedronQuadrature(Rule rule, int interpolation_order)
{
    switch (rule)
    {
        case Rule::OnePoint:
        {
            w = {1.0};
            clist = {{0, 0.25, 0.25, 0.25}};
            break;
        }
        case Rule::FourPoint:
        {
            w = {0.25, 0.25, 0.25, 0.25};

            clist = {{0, 0.5854102, 0.1381966, 0.1381966},
                     {1, 0.1381966, 0.5854102, 0.1381966},
                     {2, 0.1381966, 0.1381966, 0.5854102},
                     {3, 0.1381966, 0.1381966, 0.1381966}};
            break;
        }
        case Rule::FivePoint:
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
    w |= ranges::action::transform([](auto const& wl) { return wl / 6.0; });
}
}
