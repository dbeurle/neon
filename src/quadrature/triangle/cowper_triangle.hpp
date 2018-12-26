
#pragma once

#include "quadrature/numerical_quadrature.hpp"

namespace neon
{
class cowper_triangle : public surface_quadrature
{
public:
    explicit cowper_triangle(int const minimum_degree);
};
}
