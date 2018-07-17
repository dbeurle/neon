
#include "hexahedron_quadrature.hpp"

namespace neon
{
hexahedron_quadrature::hexahedron_quadrature(point const p)
{
    switch (p)
    {
        case point::one:
        {
            m_weights = {8.0};
            m_coordinates = {{0, 0.0, 0.0, 0.0}};
            break;
        }
        case point::six:
        {
            m_weights.resize(6, 4.0 / 3.0);
            m_coordinates = {{0, 1.0, 0.0, 0.0},
                             {1, -1.0, 0.0, 0.0},
                             {2, 0.0, 1.0, 0.0},
                             {3, 0.0, -1.0, 0.0},
                             {4, 0.0, 0.0, 1.0},
                             {5, 0.0, 0.0, -1.0}};
            break;
        }
        case point::eight:
        {
            m_weights.resize(8, 1.0);
            auto const qp = 1.0 / std::sqrt(3.0);
            m_coordinates = {{0, -qp, -qp, -qp},
                             {1, qp, -qp, -qp},
                             {2, qp, qp, -qp},
                             {3, -qp, qp, -qp},
                             {4, -qp, -qp, qp},
                             {5, qp, -qp, qp},
                             {6, qp, qp, qp},
                             {7, -qp, qp, qp}};
            break;
        }
        case point::twentyseven:
        {
            m_weights = {125.0 / 729.0, 125.0 / 729.0, 125.0 / 729.0, 125.0 / 729.0, 200.0 / 729.0,
                         200.0 / 729.0, 200.0 / 729.0, 200.0 / 729.0, 320.0 / 729.0, 200.0 / 729.0,
                         200.0 / 729.0, 200.0 / 729.0, 200.0 / 729.0, 320.0 / 729.0, 320.0 / 729.0,
                         320.0 / 729.0, 320.0 / 729.0, 512.0 / 729.0, 125.0 / 729.0, 125.0 / 729.0,
                         125.0 / 729.0, 125.0 / 729.0, 200.0 / 729.0, 200.0 / 729.0, 200.0 / 729.0,
                         200.0 / 729.0, 320.0 / 729.0};

            auto const alpha = std::sqrt(3.0 / 5.0);

            m_coordinates =
                {{0, -alpha, -alpha, -alpha}, {1, alpha, -alpha, -alpha}, {2, alpha, alpha, -alpha},
                 {3, -alpha, alpha, -alpha},  {4, 0.0, -alpha, -alpha},   {5, alpha, 0.0, -alpha},
                 {6, 0.0, alpha, -alpha},     {7, -alpha, 0.0, -alpha},   {8, 0.0, 0.0, -alpha},
                 {9, -alpha, -alpha, 0.0},    {10, alpha, -alpha, 0.0},   {11, alpha, alpha, 0.0},
                 {12, -alpha, alpha, 0.0},    {13, 0.0, -alpha, 0.0},     {14, alpha, 0.0, 0.0},
                 {15, 0.0, alpha, 0.0},       {16, -alpha, 0.0, 0.0},     {17, 0.0, 0, 0.0},
                 {18, -alpha, -alpha, alpha}, {19, alpha, -alpha, alpha}, {20, alpha, alpha, alpha},
                 {21, -alpha, alpha, alpha},  {22, 0.0, -alpha, alpha},   {23, alpha, 0.0, alpha},
                 {24, 0.0, alpha, alpha},     {25, -alpha, 0.0, alpha},   {26, 0.0, 0.0, alpha}};
            break;
        }
    }
}
}
