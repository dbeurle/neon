
#include "triangle_quadrature.hpp"

#include <algorithm>

namespace neon
{
triangle_quadrature::triangle_quadrature(point const p)
{
    switch (p)
    {
        case point::one:
        {
            m_weights = {1.0};
            m_coordinates = {{0, 1.0 / 3.0, 1.0 / 3.0}};
            break;
        }
        case point::three:
        {
            m_weights = {1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0};
            m_coordinates = {{0, 2.0 / 3.0, 1.0 / 6.0},
                             {1, 1.0 / 6.0, 2.0 / 3.0},
                             {2, 1.0 / 6.0, 1.0 / 6.0}};
            break;
        }
        case point::four:
        {
            m_weights = {-0.5625,
                         0.5208333333333333333333,
                         0.5208333333333333333333,
                         0.5208333333333333333333};
            m_coordinates = {{0, 1.0 / 3.0, 1.0 / 3.0}, {1, 0.6, 0.2}, {2, 0.2, 0.6}, {3, 0.2, 0.2}};
            break;
        }
    }
    std::transform(begin(m_weights), end(m_weights), begin(m_weights), [](auto const weight) {
        return 0.5 * weight;
    });
}
}
