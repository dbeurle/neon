
#pragma once

#include "quadrature/numerical_quadrature.hpp"

namespace neon
{
class hexahedron_quadrature : public volume_quadrature
{
public:
    explicit hexahedron_quadrature(int const minimum_degree) noexcept(false);
};
}
