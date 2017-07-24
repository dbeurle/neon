
#include "HexahedronQuadrature.hpp"

namespace neon
{
HexahedronQuadrature::HexahedronQuadrature(Rule rule, int interpolationOrder)
{
    switch (rule)
    {
        case Rule::OnePoint:
        {
            w = {8.0};
            clist = {{0, 0.0, 0.0, 0.0}};
            break;
        }
        case Rule::EightPoint:
        {
            w.resize(8, 1.0);
            auto const qp = 1.0 / std::sqrt(3.0);
            clist = {{0, -qp, -qp, -qp},
                     {1, qp, -qp, -qp},
                     {2, qp, qp, -qp},
                     {3, -qp, qp, -qp},
                     {4, -qp, -qp, qp},
                     {5, qp, -qp, qp},
                     {6, qp, qp, qp},
                     {7, -qp, qp, qp}};
            break;
        }
    }
}
}
