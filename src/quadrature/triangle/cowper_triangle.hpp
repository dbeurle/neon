
#pragma once

#include "quadrature/triangle/triangle_quadrature.hpp"

namespace neon
{
class cowper_triangle : public triangle_quadrature
{
public:
    explicit cowper_triangle(int const minimum_degree);
};
}
