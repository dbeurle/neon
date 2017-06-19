
#pragma once

#include "interpolations/ShapeFunction.hpp"

#include "quadrature/NumericalQuadrature.hpp"

namespace neon
{
/**
 * SurfaceQuadrature provides outward normal computation
 * and projection methods of three dimension coordinates to a planar surface
 */
class SurfaceInterpolation : public ShapeFunction<SurfaceQuadrature>
{
public:
    SurfaceInterpolation(std::unique_ptr<SurfaceQuadrature>&& quadratureImpl)
        : ShapeFunction<SurfaceQuadrature>(std::move(quadratureImpl))
    {
    }

    Vector3 unit_outward_normal(Matrix const& nodal_coordinates) const;

    Matrix project_to_plane(Matrix const& nodal_coordinates) const;
};
}
