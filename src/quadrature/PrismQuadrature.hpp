
#pragma once

#include "NumericalQuadrature.hpp"

namespace neon
{
class PrismQuadrature : public VolumeQuadrature
{
public:
    /** Available quadrature rules for this element type */
    enum class Rule { OnePoint, SixPoint };

public:
    PrismQuadrature(Rule rule, int interpolation_order = 1);
};
}
