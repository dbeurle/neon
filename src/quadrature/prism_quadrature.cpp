
#include "prism_quadrature.hpp"

namespace neon
{
prism_quadrature::prism_quadrature(Rule rule)
{
    switch (rule)
    {
        case Rule::OnePoint:
        {
            // One single point in the middle of the triangle and linear line
            // interpolation also in the middle
            w = {2.0};
            clist = {{0, -1.0 / 3.0, -1.0 / 3.0, 0.0}};
            break;
        }
        case Rule::SixPoint:
        {
            // A three point stencil on each triangle and two quadrature points
            // along the line element
            w.resize(6, 1.0 / 3.0 * 1.0);

            clist = {{0, 2.0 / 3.0, 1.0 / 6.0, -1.0 / std::sqrt(3.0)},
                     {1, 1.0 / 6.0, 2.0 / 3.0, -1.0 / std::sqrt(3.0)},
                     {2, 1.0 / 6.0, 1.0 / 6.0, -1.0 / std::sqrt(3.0)},
                     {3, 2.0 / 3.0, 1.0 / 6.0, 1.0 / std::sqrt(3.0)},
                     {4, 1.0 / 6.0, 2.0 / 3.0, 1.0 / std::sqrt(3.0)},
                     {5, 1.0 / 6.0, 1.0 / 6.0, 1.0 / std::sqrt(3.0)}};

            break;
        }
    }
}
}
