
#pragma once

/// @file

#include "quadrature/numerical_quadrature.hpp"

namespace neon
{
class quadrilateral_quadrature : public surface_quadrature
{
public:
    explicit quadrilateral_quadrature(int const minimum_degree);
};
}
