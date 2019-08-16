
#include "gauss_legendre.hpp"

#include <cmath>

namespace neon::quadrature::quadrilateral
{
gauss_legendre::gauss_legendre(int const minimum_degree)
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
            m_degree = 3;
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
            m_degree = 5;
            m_weights = {25.0 / 81.0,
                         25.0 / 81.0,
                         25.0 / 81.0,
                         25.0 / 81.0,
                         40.0 / 81.0,
                         40.0 / 81.0,
                         40.0 / 81.0,
                         40.0 / 81.0,
                         64.0 / 81.0};
            m_coordinates = {{0, -std::sqrt(3.0 / 5.0), -std::sqrt(3.0 / 5.0)},
                             {1, std::sqrt(3.0 / 5.0), -std::sqrt(3.0 / 5.0)},
                             {2, std::sqrt(3.0 / 5.0), std::sqrt(3.0 / 5.0)},
                             {3, -std::sqrt(3.0 / 5.0), std::sqrt(3.0 / 5.0)},
                             {4, 0.0, -std::sqrt(3.0 / 5.0)},
                             {5, std::sqrt(3.0 / 5.0), 0.0},
                             {6, 0.0, std::sqrt(3.0 / 5.0)},
                             {7, -std::sqrt(3.0 / 5.0), 0.0},
                             {8, 0.0, 0.0}};
            break;
        }
        default:
        {
            throw std::domain_error("quadrature scheme greater than degree 5 not yet available");
        }
    }
}
}
