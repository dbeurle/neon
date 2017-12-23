
#pragma once

#include "ShapeFunction.hpp"
#include "quadrature/PrismQuadrature.hpp"

namespace neon
{
/**
 * prism6 is responsible for computing the shape functions for an prism or
 * wedge shaped six noded element.
 * The shape functions and nodal orderings are taken from \cite{Hughes2012}
 */
class prism6 : public VolumeInterpolation
{
public:
    explicit prism6(PrismQuadrature::Rule rule);

    virtual int nodes() const override final { return 6; }

    double compute_measure(matrix3x const& nodal_coordinates) const;

protected:
    /**
     * Initialise the shape functions to the following polynomials
     * \f{align*}{
     * N_1(r, s, \xi) &= \frac{1}{2}(1-\xi)r \\
     * N_2(r, s, \xi) &= \frac{1}{2}(1+\xi)s \\
     * N_3(r, s, \xi) &= \frac{1}{2}(1+\xi)t \\
     * N_4(r, s, \xi) &= \frac{1}{2}(1-\xi)r \\
     * N_5(r, s, \xi) &= \frac{1}{2}(1-\xi)s \\
     * N_6(r, s, \xi) &= \frac{1}{2}(1+\xi)r \\
     * \f}
     */
    void precompute_shape_functions();
};
}
