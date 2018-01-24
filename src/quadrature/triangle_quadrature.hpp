
#pragma once

#include "numerical_quadrature.hpp"

namespace neon
{
class triangle_quadrature : public surface_quadrature
{
public:
    /** Available quadrature rules for this element type */
    enum class Rule { OnePoint, ThreePoint, FourPoint };

public:
    explicit triangle_quadrature(Rule rule);
};
}
