
#pragma once

#include "ShapeFunction.hpp"
#include "quadrature/TetrahedronQuadrature.hpp"

namespace neon
{
/**
 * Tetrahedron4 is an isoparametric tetrahedral 4 node element with analytical
 * integration.
 */
class Tetrahedron4 : public VolumeInterpolation
{
public:
    Tetrahedron4(TetrahedronQuadrature::Rule rule = TetrahedronQuadrature::Rule::OnePoint);

    virtual int nodes() const override final { return 4; }

protected:
    /**
     * Initialize the shape functions to the following polynomials
     * \f{align*}{
     * N_1(r, s, t) &= r \\
     * N_2(r, s, t) &= s \\
     * N_3(r, s, t) &= t \\
     * N_4(r, s, t) &= 1 - r - s - t
     * \f}
     */
    void precompute_shape_functions();
};
}
