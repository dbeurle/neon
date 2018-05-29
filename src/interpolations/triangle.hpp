
#pragma once

#include "shape_function.hpp"
#include "quadrature/triangle_quadrature.hpp"

namespace neon
{
/** Triangular 3 node element with analytic integration */
class triangle3 : public surface_interpolation
{
public:
    triangle3(triangle_quadrature::point const p);

    virtual int nodes() const override final { return 3; }

    double compute_measure(matrix const& nodal_coordinates) const;

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

class triangle6 : public surface_interpolation
{
public:
    triangle6(triangle_quadrature::point const p);

    virtual int nodes() const override final { return 6; }

    /**
     * Compute the area of the element using Gaussian integration
     * @param nodal_coordinates element nodal coordinates
     * @return element face area
     */
    double compute_measure(matrix const& nodal_coordinates);

protected:
    void precompute_shape_functions();
};
}
