
#pragma once

#include "NumericalQuadrature.hpp"

namespace neon
{
class HexahedronQuadrature : public VolumeQuadrature
{
public:
    /** Available quadrature rules for this element type */
    enum class Rule { OnePoint, EightPoint };

public:
    HexahedronQuadrature(Rule rule, int interpolationOrder = 1);
};
}
