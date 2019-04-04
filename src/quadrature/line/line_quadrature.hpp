
#pragma once

/// @file

#include "quadrature/numerical_quadrature.hpp"

namespace neon
{
class line_quadrature : public numerical_quadrature<double>
{
public:
    /// Available quadrature rules
    enum class point {
        /// First order
        one,
        /// Third order
        two,
        /// Fifth order
        three,
        /// Seventh order
        four
    };

public:
    explicit line_quadrature(point const p);
};
}
