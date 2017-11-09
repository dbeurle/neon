
#include "PrismQuadrature.hpp"

namespace neon
{
PrismQuadrature::PrismQuadrature(Rule rule)
{
    switch (rule)
    {
        case Rule::OnePoint:
        {
            w = {4.0};

            clist = {{0, -1.0 / 3.0, -1.0 / 3.0, 0.0}};

            break;
        }
        case Rule::SixPoint:
        {
            w.resize(6, 1.0);

            clist = {{0, 2.0 / 3.0, 1.0 / 6.0, 1.0 / 6.0},
                     {1, 2.0 / 3.0, 1.0 / 6.0, 1.0 / 6.0},
                     {2, 1.0 / 6.0, 2.0 / 3.0, 1.0 / 6.0},
                     {3, 1.0 / 6.0, 2.0 / 3.0, 1.0 / 6.0},
                     {4, -1.0 / std::sqrt(3.0), -1.0 / std::sqrt(3.0), -1.0 / std::sqrt(3.0)},
                     {5, 1.0 / std::sqrt(3.0), 1.0 / std::sqrt(3.0), 1.0 / std::sqrt(3.0)}};

            break;

            // TODO Check if this quadrature scheme is useful
            // m_quadraturePoints = 3;
            // m_quadratureWeight = 1.0;
            // m_quadratureCoordinates.resize(3, 4);
            // coordinates << 0.483163247594393, -0.605498860309242, -0.605498860309242, //
            //		-0.741581623797196, 0.469416096821288,
            //-0.530583903178712,
            //
            //		0.0, 0.0,
        }
    }
}
}
