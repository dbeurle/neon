
#pragma once

#include "NumericalQuadrature.hpp"

namespace neon
{
class HexahedronQuadrature : public VolumeQuadrature
{
public:
    /** Available quadrature rules for this element type */
    enum class Rule {
        OnePoint,
        SixPoint, // Fourth order accuracy
        EightPoint
    };

public:
    HexahedronQuadrature(Rule rule);
};
}
