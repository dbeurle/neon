
#pragma once

#include "quadrature/numerical_quadrature.hpp"

namespace neon::quadrature::tetrahedron
{
class jinyun : public volume_quadrature
{
public:
    explicit jinyun(int const minimum_degree) noexcept(false);
};
}
