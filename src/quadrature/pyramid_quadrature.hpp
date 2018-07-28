
#pragma once

#include "numerical_quadrature.hpp"

namespace neon
{
class pyramid_quadrature : public volume_quadrature
{
public:
    /// Available quadrature rules for this element type
    enum class point {
        /// First order
        one,
        /// Third order
        eight,
        /// Fifth order
        twenty_seven,
        /// Seventh order
        sixty_four
    };

public:
    explicit pyramid_quadrature(point const p);
};
}
