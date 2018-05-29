
#include "line_quadrature.hpp"

namespace neon
{
line_quadrature::line_quadrature(point const p)
{
    switch (p)
    {
        case point::one:
        {
            w = {2.0};
            clist = {{0, 0.0}};
            break;
        }
        case point::two:
        {
            w = {1.0, 1.0};

            clist = {{0, -1.0 / std::sqrt(3.0)}, {1, 1.0 / std::sqrt(3.0)}};
            break;
        }
        case point::three:
        {
            w = {5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0};

            clist = {{0, -std::sqrt(3.0 / 5.0)}, {1, 0.0}, {2, std::sqrt(3.0 / 5.0)}};
            break;
        }
    }
}
}
