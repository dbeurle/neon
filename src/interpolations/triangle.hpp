
#pragma once

/// @file

#include "shape_function.hpp"
#include "quadrature/triangle/triangle_quadrature.hpp"

namespace neon
{
/// Triangular 3 node element with analytic integration
class triangle3 : public surface_interpolation
{
public:
    triangle3(triangle_quadrature::point const p);

    double compute_measure(matrix const& nodal_coordinates) const;

protected:
    /// Initialize the shape functions
    void precompute_shape_functions();
};

/// A six node triangle with midside nodes and quadrature shape functions.
/// \sa Hughes2012
class triangle6 : public surface_interpolation
{
public:
    triangle6(triangle_quadrature::point const p);

    /// Compute the area of the element using numerical integration
    /// \param nodal_coordinates element nodal coordinates
    /// \return surface area
    double compute_measure(matrix const& nodal_coordinates);

protected:
    void precompute_shape_functions();
};
}
