
#pragma once

#include "ShapeFunction.hpp"
#include "quadrature/TriangleQuadrature.hpp"

namespace neon
{
/** Triangular 3 node element with analytic integration */
class Triangle3 : public SurfaceInterpolation
{
public:
    Triangle3(TriangleQuadrature::Rule rule);

    virtual int nodes() const override final { return 3; }

    double compute_measure(Matrix const& nodal_coordinates) const;

protected:
    /**
     * Initialize the shape functions to the following polynomials
     * \f{align*}{
     * N_1(r, s, t) &= r \\
     * N_2(r, s, t) &= s \\
     * N_3(r, s, t) &= 1 - r - s
     * \f}
     */
    void precompute_shape_functions();
};
}
