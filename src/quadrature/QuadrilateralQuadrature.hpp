
#pragma once

#include "NumericalQuadrature.hpp"

namespace neon
{
class QuadrilateralQuadrature : public SurfaceQuadrature
{
public:
    /** Available quadrature rules for this element type */
    enum class Rule { OnePoint, FourPoint, NinePoint };

public:
    explicit QuadrilateralQuadrature(Rule rule);
};
}
