
#include "gauss_legendre_line.hpp"

#include <cmath>

namespace neon
{
gauss_legendre_line::gauss_legendre_line(int const minimum_degree)
{
    switch (minimum_degree)
    {
        case 1:
        {
            m_degree = 1;
            m_weights = {2.0};
            m_coordinates = {{0, 0.0}};
            break;
        }
        case 2:
        case 3:
        {
            m_degree = 3;
            m_weights = {1.0, 1.0};
            m_coordinates = {{0, -1.0 / std::sqrt(3.0)}, {1, 1.0 / std::sqrt(3.0)}};
            break;
        }
        case 4:
        case 5:
        {
            m_degree = 5;
            m_weights = {5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0};
            m_coordinates = {{0, -std::sqrt(3.0 / 5.0)}, {1, 0.0}, {2, std::sqrt(3.0 / 5.0)}};
            break;
        }
        case 6:
        case 7:
        {
            m_degree = 7;
            m_weights = {0.347854845137454, 0.652145154862546, 0.652145154862546, 0.347854845137454};
            m_coordinates = {{0, -0.861136311594053},
                             {1, -0.339981043584856},
                             {2, 0.339981043584856},
                             {3, 0.861136311594053}};
            break;
        }
        default:
        {
            throw std::domain_error("Gauss-Legendre quadrature currently supports up to 7th "
                                    "degree");
        }
    }
}
}
