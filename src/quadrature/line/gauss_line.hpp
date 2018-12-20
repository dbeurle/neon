
#pragma once

#include "quadrature/numerical_quadrature.hpp"

namespace neon
{
class gauss_line : public numerical_quadrature<double>
{
public:
    explicit gauss_line(int const minimum_degree);
};
}
