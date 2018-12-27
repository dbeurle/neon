
#pragma once

#include "shape_function.hpp"

namespace neon
{
/// Line element with two nodes with Gaussian integration
class line2 : public line_interpolation
{
public:
    explicit line2();

    /// Evaluate the shape functions at the natural coordinate
    virtual auto evaluate(coordinate_type const& coordinate) const noexcept(false)
        -> value_type override final;
};

/// Isoparameteric line element with three nodes using Gaussian quadrature
/**
 * Implements shape functions according to the following polynomials
 * \f{align*}{
 * N_1(r, s, t) &= \frac{1}{2} \xi \left(\xi - 1 \right) \\
 * N_2(r, s, t) &= 1 - \xi^2 \\
 * N_3(r, s, t) &= \frac{1}{2} \xi \left(\xi + 1 \right)
 * \f}
 */
class line3 : public line_interpolation
{
public:
    explicit line3();

    /// Evaluate the shape functions at the natural coordinate
    virtual auto evaluate(coordinate_type const& coordinate) const noexcept(false)
        -> value_type override final;
};
}
