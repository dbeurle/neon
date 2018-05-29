
#pragma once

#include "numerical_quadrature.hpp"

namespace neon
{
class quadrilateral_quadrature : public surface_quadrature
{
public:
    /// Available quadrature rules for this element type
    enum class point { one, four, nine };

public:
    explicit quadrilateral_quadrature(point const p);
};
}
