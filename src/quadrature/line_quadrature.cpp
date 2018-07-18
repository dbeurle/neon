
#include "line_quadrature.hpp"

namespace neon
{
line_quadrature::line_quadrature(point const p)
{
    switch (p)
    {
        case point::one:
        {
            m_weights = {2.0};
            m_coordinates = {{0, 0.0}};
            break;
        }
        case point::two:
        {
            m_weights = {1.0, 1.0};

            m_coordinates = {{0, -1.0 / std::sqrt(3.0)}, {1, 1.0 / std::sqrt(3.0)}};
            break;
        }
        case point::three:
        {
            m_weights = {5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0};

            m_coordinates = {{0, -std::sqrt(3.0 / 5.0)}, {1, 0.0}, {2, std::sqrt(3.0 / 5.0)}};
            break;
        }
        case point::four:
        {
            m_weights = {0.347854845137454, 0.652145154862546, 0.652145154862546, 0.347854845137454};
            m_coordinates = {{0, -0.861136311594053},
                             {1, -0.339981043584856},
                             {2, 0.339981043584856},
                             {3, 0.861136311594053}};
            break;
        }
    }
}
}
