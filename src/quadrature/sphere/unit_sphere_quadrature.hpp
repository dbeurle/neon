
#pragma once

/// @file

#include "quadrature/numerical_quadrature.hpp"

namespace neon
{
/// unit_sphere_quadrature implements the symmetric quadrature rules
class unit_sphere_quadrature : public volume_quadrature
{
public:
    /// Quadrature nomenclature from @cite Ehret2010
    enum class point { BO21, BO33, BO61, FM900 };

public:
    /// Fill the quadrature coordinates and weightings
    unit_sphere_quadrature(point const p);

protected:
    void precompute_coordinates();
};
}
