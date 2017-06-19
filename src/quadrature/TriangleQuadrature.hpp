
#pragma once

#include "NumericalQuadrature.hpp"

namespace neon
{
class TriangleQuadrature : public SurfaceQuadrature
{
public:
    /** Available quadrature rules for this element type */
    enum class Rule { OnePoint, ThreePoint, FourPoint };

public:
    TriangleQuadrature(Rule rule, int interpolation_order = 1);
};
}
