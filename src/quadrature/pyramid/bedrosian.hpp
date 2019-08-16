
#pragma once

/// @file

#include "quadrature/numerical_quadrature.hpp"

namespace neon::quadrature::pyramid
{
/// Pyramid quadrature scheme contructed by combining quadrilateral
/// based on the schemes by Hammer et.al @cite Bedrosian1992
class bedrosian : public volume_quadrature
{
public:
    explicit bedrosian(int const minimum_degree);
};
}
