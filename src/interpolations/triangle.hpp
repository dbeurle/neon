
#pragma once

#include "shape_function.hpp"
#include "quadrature/triangle_quadrature.hpp"

namespace neon
{
/// Triangle three node element with constant derivatives with shape functions
/// given by
/// \f{align*}{
/// N_1(r, s, t) &= r \\
/// N_2(r, s, t) &= s \\
/// N_3(r, s, t) &= 1 - r - s
/// \f}
/// \sa Hughes2012
class triangle3 : public surface_interpolation
{
public:
    triangle3(triangle_quadrature::point const p);

    int nodes() const override final { return 3; }

    double compute_measure(matrix const& nodal_coordinates) const;

protected:
    /**
     * Initialize the shape functions to the following polynomials
     */
    void precompute_shape_functions();
};

/// A six node triangle with midside nodes and quadrature shape functions.
/// \sa Hughes2012
class triangle6 : public surface_interpolation
{
public:
    triangle6(triangle_quadrature::point const p);

    int nodes() const override final { return 6; }

    /// Compute the area of the element using numerical integration
    /// \param nodal_coordinates element nodal coordinates
    /// \return surface area
    double compute_measure(matrix const& nodal_coordinates);

protected:
    void precompute_shape_functions();
};
}
