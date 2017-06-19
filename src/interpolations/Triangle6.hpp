
#pragma once

#include "SurfaceInterpolation.hpp"
#include "quadrature/TriangleQuadrature.hpp"

namespace neon
{
class Triangle6 : public SurfaceInterpolation
{
public:
    Triangle6(TriangleQuadrature::Rule rule);

    virtual int nodes() const override final { return 6; }

    /**
     * Compute the area using Gaussian integration.
     * @param nodal coordinates
     * @return Element face area
     */
    double compute_measure(const Matrix& nodal_coordinates);

protected:
    void precompute_shape_functions();
};
}
