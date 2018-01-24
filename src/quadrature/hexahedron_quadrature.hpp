
#pragma once

#include "numerical_quadrature.hpp"

namespace neon
{
class hexahedron_quadrature : public volume_quadrature
{
public:
    /** Available quadrature rules for this element type */
    enum class Rule {
        OnePoint,
        SixPoint,        // Fourth order accuracy
        EightPoint,      // Fourth order accuracy
        TwentySevenPoint // Six order accuracy
    };

public:
    hexahedron_quadrature(Rule const rule);
};
}
