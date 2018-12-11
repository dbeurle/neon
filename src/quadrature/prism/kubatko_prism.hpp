
#pragma once

#include "quadrature/numerical_quadrature.hpp"

namespace neon
{
/// Prism schemes for points that are inside the prism with positive weights
/// \cite Kubatko2013
class kubatko_prism : public volume_quadrature
{
public:
    explicit kubatko_prism(int const minimum_degree);
};
}
