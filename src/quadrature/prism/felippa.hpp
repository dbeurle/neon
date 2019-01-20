
#pragma once

/// @file

#include "quadrature/numerical_quadrature.hpp"

namespace neon::quadrature::prism
{
/// Prism schemes for points that are inside the prism with positive weights
/// @cite Kubatko2013
class felippa : public volume_quadrature
{
public:
    /// Construct using a minimum degree and will throw a std::domain_error
    /// if the minimum degree cannot be met
    explicit felippa(int const minimum_degree) noexcept(false);
};
}
