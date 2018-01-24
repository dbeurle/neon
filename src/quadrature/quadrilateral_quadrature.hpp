
#pragma once

#include "numerical_quadrature.hpp"

namespace neon
{
class quadrilateral_quadrature : public surface_quadrature
{
public:
    /** Available quadrature rules for this element type */
    enum class Rule { OnePoint, FourPoint, NinePoint };

public:
    explicit quadrilateral_quadrature(Rule rule);
};
}
