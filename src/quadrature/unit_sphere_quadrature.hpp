
#pragma once

#include "numerical_quadrature.hpp"

namespace neon
{
/** unit_sphere_quadrature implements the symmetric quadrature rules */
class unit_sphere_quadrature : public volume_quadrature
{
public:
    /** Quadrature nomenclature from \cite Ehret2010 */
    enum class Rule { BO21, BO33, FM900 };

public:
    /** Fill the quadrature coordinates and weightings */
    unit_sphere_quadrature(Rule const rule);

protected:
    void precompute_coordinates();
};
}
