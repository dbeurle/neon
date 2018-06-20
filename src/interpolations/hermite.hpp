
#pragma once

#include "shape_function.hpp"
#include "quadrature/line_quadrature.hpp"

namespace neon
{
class hermite : public line_interpolation
{
public:
    hermite(line_quadrature::point const p);

    virtual int nodes() const override final { return 2; }

    /// \return element length
    double compute_measure(matrix const& nodal_coordinates) const;

protected:
    /**
     * Initialize the shape functions to the following polynomials
     * \f{align*}{
     * N_1(\xi) &= \frac{1}{2}\left(1 - \xi \right) \\
     * N_2(\xi) &= \frac{1}{2}\left(1 + \xi \right) \\
     * N_3(\xi) &= \frac{1}{2}\left(1 + \xi \right) \\
     * N_4(\xi) &= \frac{1}{2}\left(1 + \xi \right)
     * \f}
     */
    void precompute_shape_functions();
};
}
