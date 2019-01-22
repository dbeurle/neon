
#pragma once

/// @file

#include "quadrature/numerical_quadrature.hpp"

namespace neon::quadrature::prism
{
/// Prism schemes for points that are inside the prism with positive weights
/// @cite Kubatko2013
class kubatko : public volume_quadrature
{
public:
    /// Construct using a minimum degree and will throw a std::domain_error
    /// if the scheme is not available or is not able to be met
    explicit kubatko(int const minimum_degree) noexcept(false);
};
}
