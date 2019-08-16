
#pragma once

#include "quadrature/numerical_quadrature.hpp"

namespace neon::quadrature::tetrahedron
{
class witherden_vincent : public volume_quadrature
{
public:
    explicit witherden_vincent(int const minimum_degree);
};
}
