
#pragma once

#include "numerical_quadrature.hpp"

namespace neon
{
class prism_quadrature : public volume_quadrature
{
public:
    /** Available quadrature rules for this element type */
    enum class Rule { OnePoint, SixPoint };

public:
    explicit prism_quadrature(Rule rule);
};
}
