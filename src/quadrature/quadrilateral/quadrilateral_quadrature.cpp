
#include "quadrilateral_quadrature.hpp"

namespace neon
{
quadrilateral_quadrature::quadrilateral_quadrature(point const p)
{
    switch (p)
    {
        case point::one:
        {
            m_weights = {4.0};
            m_coordinates = {{0, 0.0, 0.0}};
            break;
        }
        case point::four:
        {
            m_weights = {1.0, 1.0, 1.0, 1.0};
            m_coordinates = {{0, -1.0 / std::sqrt(3.0), -1.0 / std::sqrt(3.0)},
                             {1, 1.0 / std::sqrt(3.0), -1.0 / std::sqrt(3.0)},
                             {2, 1.0 / std::sqrt(3.0), 1.0 / std::sqrt(3.0)},
                             {3, -1.0 / std::sqrt(3.0), 1.0 / std::sqrt(3.0)}};
            break;
        }
        case point::nine:
        {
            m_weights = {25.0 / 81.0,
                         25.0 / 81.0,
                         25.0 / 81.0,
                         25.0 / 81.0,
                         40.0 / 81.0,
                         40.0 / 81.0,
                         40.0 / 81.0,
                         40.0 / 81.0,
                         64.0 / 81.0};

            auto const coordinate = std::sqrt(3.0 / 5.0);

            m_coordinates = {{0, -coordinate, -coordinate},
                             {1, coordinate, -coordinate},
                             {2, coordinate, coordinate},
                             {3, -coordinate, coordinate},
                             {4, 0.0, -coordinate},
                             {5, coordinate, 0.0},
                             {6, 0.0, coordinate},
                             {7, -coordinate, 0.0},
                             {8, 0.0, 0.0}};
            break;
        }
    }
}
}
