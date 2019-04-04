
#pragma once

/// @file

#include "quadrature/numerical_quadrature.hpp"

namespace neon
{
class hexahedron_quadrature : public volume_quadrature
{
public:
    /// Available quadrature rules for this element type
    enum class point {
        /// Centroid first order
        one,
        /// Fourth order accuracy
        six,
        /// Fourth order accuracy (recommended over six)
        eight,
        /// Sixth order accuracy
        twentyseven
    };

public:
    hexahedron_quadrature(point const p);
};
}
