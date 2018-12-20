
#pragma once

#include "quadrature/numerical_quadrature.hpp"

namespace neon
{
class tetrahedron_jinyun : public volume_quadrature
{
public:
    explicit tetrahedron_jinyun(int const minimum_degree);
};
}
