
#pragma once

#include "ShapeFunction.hpp"
#include "quadrature/LineQuadrature.hpp"

namespace neon
{
/** Line 2 node element with Gaussian integration */
class Line2 : public LineInterpolation
{
public:
    Line2(LineQuadrature::Rule rule);

    virtual int nodes() const override final { return 2; }

    double compute_measure(matrix const& nodal_coordinates) const;

protected:
    /**
     * Initialize the shape functions to the following polynomials
     * \f{align*}{
     * N_1(\xi) &= \frac{1}{2}\left(1 - \xi \right) \\
     * N_2(\xi) &= \frac{1}{2}\left(1 + \xi \right)
     * \f}
     */
    void precompute_shape_functions();
};

/** Line 3 node element with Gaussian integration */
class Line3 : public LineInterpolation
{
public:
    Line3(LineQuadrature::Rule rule);

    virtual int nodes() const override final { return 3; }

    double compute_measure(matrix const& nodal_coordinates) const;

protected:
    /**
     * Initialize the shape functions to the following polynomials
     * \f{align*}{
     * N_1(\xi) &= \frac{1}{2} \xi \left(\xi - 1 \right) \\
     * N_2(\xi) &= 1 - \xi^2 \\
     * N_3(\xi) &= \frac{1}{2} \xi \left(\xi + 1 \right)
     * \f}
     */
    void precompute_shape_functions();
};
}
