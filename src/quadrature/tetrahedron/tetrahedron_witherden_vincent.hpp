
#pragma once

#include "quadrature/numerical_quadrature.hpp"

namespace neon
{
class tetrahedron_witherden_vincent : public volume_quadrature
{
    explicit tetrahedron_witherden_vincent(int const minimum_degree);
};
}
