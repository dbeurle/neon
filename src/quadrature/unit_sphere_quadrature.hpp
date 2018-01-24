
#pragma once

#include "numerical_quadrature.hpp"

namespace neon
{
/** unit_sphere_quadrature implements the symmetric quadrature rules */
class unit_sphere_quadrature : public volume_quadrature
{
public:
    /**
     * Quadrature nomenclature from
     * Numerical integration on the sphere and its effect on the material
     * symmetry of constitutive equations - A comparative study
     * by Ehret et.al.
     */
    enum class Rule { BO21, BO33, FM900 };

public:
    /** Fill the quadrature coordinates and weightings */
    unit_sphere_quadrature(Rule const rule);

protected:
    void precompute_coordinates();
};
}
