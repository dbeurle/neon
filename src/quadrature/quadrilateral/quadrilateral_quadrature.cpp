
#include "quadrilateral_quadrature.hpp"

#include <cmath>

namespace neon
{
quadrilateral_quadrature::quadrilateral_quadrature(int const minimum_degree)
{
    switch (minimum_degree)
    {
        case 1:
        {
            m_degree = 1;
            m_weights = {4.0};
            m_coordinates = {{0, 0.0, 0.0}};
            break;
        }
        case 2:
        case 3:
        {
            m_weights = {1.0, 1.0, 1.0, 1.0};

            m_coordinates = {{0, -1.0 / std::sqrt(3.0), -1.0 / std::sqrt(3.0)},
                             {1, 1.0 / std::sqrt(3.0), -1.0 / std::sqrt(3.0)},
                             {2, 1.0 / std::sqrt(3.0), 1.0 / std::sqrt(3.0)},
                             {3, -1.0 / std::sqrt(3.0), 1.0 / std::sqrt(3.0)}};
            break;
        }
        case 4:
        case 5:
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

            auto const a = std::sqrt(3.0 / 5.0);

            m_coordinates = {{0, -a, -a},
                             {1, a, -a},
                             {2, a, a},
                             {3, -a, a},
                             {4, 0.0, -a},
                             {5, a, 0.0},
                             {6, 0.0, a},
                             {7, -a, 0.0},
                             {8, 0.0, 0.0}};
            break;
        }
    }
}
}
