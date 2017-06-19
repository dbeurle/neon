
#pragma once

#include "SurfaceInterpolation.hpp"
#include "quadrature/TriangleQuadrature.hpp"

namespace neon
{
/** Triangular 3 node element with analytic integration */
class Triangle3 : public SurfaceInterpolation
{
public:
    Triangle3(TriangleQuadrature::Rule rule);

    virtual int nodes() const override final { return 3; }

    double compute_measure(Matrix const& nodal_coordinates);

protected:
    void precompute_shape_functions();
};
}
