
#pragma once

#include "numerical_quadrature.hpp"

namespace neon
{
class triangle_quadrature : public surface_quadrature
{
public:
    /** Available quadrature rules for this element type */
    enum class point { one, three, four };

public:
    explicit triangle_quadrature(point const p);
};
}
