
#pragma once

#include "numerical_quadrature.hpp"

namespace neon
{
class line_quadrature : public numerical_quadrature<double>
{
public:
    /** Available quadrature rules for this element type */
    enum class Rule { OnePoint, TwoPoint, ThreePoint };

public:
    explicit line_quadrature(Rule rule);
};
}
