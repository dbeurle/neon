
#pragma once

#include "quadrature/numerical_quadrature.hpp"

namespace neon
{
/// \cite Bedrosian1992
class prism_quadrature : public volume_quadrature
{
public:
    /// Available quadrature rules for this element type
    enum class point { one, six, nine };

public:
    explicit prism_quadrature(point const p);
};
}
