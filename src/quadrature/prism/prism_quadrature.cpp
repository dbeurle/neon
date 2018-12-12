
#include "prism_quadrature.hpp"

namespace neon
{
prism_quadrature::prism_quadrature(point const p)
{
    switch (p)
    {
        case point::one:
        {
            // One single point in the middle of the triangle and linear line
            // interpolation also in the middle
            m_weights = {0.5 * 2.0};
            m_coordinates = {{0, 1.0 / 3.0, 1.0 / 3.0, 0.0}};
            break;
        }
        case point::six:
        {
            // A three point stencil on each triangle and two quadrature points
            // along the line element
            m_weights.resize(6, 1.0 / 6.0 * 1.0);

            m_coordinates = {{0, 2.0 / 3.0, 1.0 / 6.0, -1.0 / std::sqrt(3.0)},
                             {1, 1.0 / 6.0, 2.0 / 3.0, -1.0 / std::sqrt(3.0)},
                             {2, 1.0 / 6.0, 1.0 / 6.0, -1.0 / std::sqrt(3.0)},
                             {3, 2.0 / 3.0, 1.0 / 6.0, 1.0 / std::sqrt(3.0)},
                             {4, 1.0 / 6.0, 2.0 / 3.0, 1.0 / std::sqrt(3.0)},
                             {5, 1.0 / 6.0, 1.0 / 6.0, 1.0 / std::sqrt(3.0)}};

            break;
        }
        case point::nine:
        {
            // Three points on each triangle and three points along the line element
            m_weights.resize(9, 1.0 / 6.0);
            m_weights.at(0) *= 5.0 / 9.0;
            m_weights.at(1) *= 5.0 / 9.0;
            m_weights.at(2) *= 5.0 / 9.0;

            m_weights.at(3) *= 8.0 / 9.0;
            m_weights.at(4) *= 8.0 / 9.0;
            m_weights.at(5) *= 8.0 / 9.0;

            m_weights.at(6) *= 5.0 / 9.0;
            m_weights.at(7) *= 5.0 / 9.0;
            m_weights.at(8) *= 5.0 / 9.0;

            m_coordinates = {{0, 2.0 / 3.0, 1.0 / 6.0, -std::sqrt(3.0 / 5.0)},
                             {1, 1.0 / 6.0, 2.0 / 3.0, -std::sqrt(3.0 / 5.0)},
                             {2, 1.0 / 6.0, 1.0 / 6.0, -std::sqrt(3.0 / 5.0)},
                             //
                             {3, 2.0 / 3.0, 1.0 / 6.0, 0.0},
                             {4, 1.0 / 6.0, 2.0 / 3.0, 0.0},
                             {5, 1.0 / 6.0, 1.0 / 6.0, 0.0},
                             //
                             {6, 2.0 / 3.0, 1.0 / 6.0, std::sqrt(3.0 / 5.0)},
                             {7, 1.0 / 6.0, 2.0 / 3.0, std::sqrt(3.0 / 5.0)},
                             {8, 1.0 / 6.0, 1.0 / 6.0, std::sqrt(3.0 / 5.0)}};

            break;
        }
    }
}
}
