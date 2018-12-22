
#pragma once

/// @file

#include "quadrature/numerical_quadrature.hpp"

namespace neon
{
/// Pyramid quadrature scheme contructed by combining quadrilateral
/// based on the schemes by Hammer et.al @cite Bedrosian1992
class bedrosian_pyramid : public volume_quadrature
{
public:
    explicit bedrosian_pyramid(int const minimum_degree);
};
}
