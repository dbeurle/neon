
#pragma once

#include "ShapeFunction.hpp"
#include "quadrature/TriangleQuadrature.hpp"

namespace neon
{
class Triangle6 : public SurfaceInterpolation
{
public:
    Triangle6(TriangleQuadrature::Rule rule);

    virtual int nodes() const override final { return 6; }

    /**
     * Compute the area of the element using Gaussian integration
     * @param nodal_coordinates element nodal coordinates
     * @return element face area
     */
    double compute_measure(Matrix const& nodal_coordinates);

protected:
    void precompute_shape_functions();
};
}
