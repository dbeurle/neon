
#pragma once

#include "NumericalQuadrature.hpp"

namespace neon
{
class LineQuadrature : public NumericalQuadrature<double>
{
public:
    /** Available quadrature rules for this element type */
    enum class Rule { OnePoint, TwoPoint, ThreePoint };

public:
    explicit LineQuadrature(Rule rule);
};
}
