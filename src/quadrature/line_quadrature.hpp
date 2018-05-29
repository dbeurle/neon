
#pragma once

#include "numerical_quadrature.hpp"

namespace neon
{
class line_quadrature : public numerical_quadrature<double>
{
public:
    /// Available quadrature rules for this element type
    enum class point { one, two, three };

public:
    explicit line_quadrature(point const p);
};
}
