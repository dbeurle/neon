
#include "QuadrilateralQuadrature.hpp"

namespace neon
{
QuadrilateralQuadrature::QuadrilateralQuadrature(Rule rule, int interpolation_order)
{
    switch (rule)
    {
        case Rule::OnePoint:
        {
            w = {4.0};
            clist = {{0, 0.0, 0.0}};
            break;
        }
        case Rule::FourPoint:
        {
            w = {1.0, 1.0, 1.0, 1.0};

            clist = {{0, -1.0 / std::sqrt(3.0), -1.0 / std::sqrt(3.0)},
                     {1, 1.0 / std::sqrt(3.0), -1.0 / std::sqrt(3.0)},
                     {2, 1.0 / std::sqrt(3.0), 1.0 / std::sqrt(3.0)},
                     {3, -1.0 / std::sqrt(3.0), 1.0 / std::sqrt(3.0)}};
            break;
        }
    }
}
}
