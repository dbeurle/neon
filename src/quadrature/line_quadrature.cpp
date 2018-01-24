
#include "line_quadrature.hpp"

namespace neon
{
line_quadrature::line_quadrature(Rule rule)
{
    switch (rule)
    {
        case Rule::OnePoint:
        {
            w = {2.0};
            clist = {{0, 0.0}};
            break;
        }
        case Rule::TwoPoint:
        {
            w = {1.0, 1.0};

            clist = {{0, -1.0 / std::sqrt(3.0)}, {1, 1.0 / std::sqrt(3.0)}};
            break;
        }
        case Rule::ThreePoint:
        {
            w = {5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0};

            clist = {{0, -std::sqrt(3.0 / 5.0)}, {1, 0.0}, {2, std::sqrt(3.0 / 5.0)}};
            break;
        }
    }
}
}
