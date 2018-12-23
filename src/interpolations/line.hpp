
#pragma once

/// @file

#include "shape_function.hpp"
#include "quadrature/line/line_quadrature.hpp"

namespace neon
{
/// Line element with two nodes with Gaussian integration
class line2 : public line_interpolation
{
public:
    explicit line2(line_quadrature::point const p);

    /// \return Element length
    double compute_measure(matrix const& nodal_coordinates) const;

protected:
    /**
     * Initialize the shape functions to the following polynomials
     * \f{align*}{
     * N_1(r, s, t) &= \frac{1}{2}\left(1 - \xi \right) \\
     * N_2(r, s, t) &= \frac{1}{2}\left(1 + \xi \right)
     * \f}
     */
    void precompute_shape_functions();
};

/// Isoparameteric line element with three nodes using Gaussian quadrature
class line3 : public line_interpolation
{
public:
    explicit line3(line_quadrature::point const p);

    /// \return Element length
    double compute_measure(matrix const& nodal_coordinates) const;

protected:
    /**
     * Initialize the shape functions to the following polynomials
     * \f{align*}{
     * N_1(r, s, t) &= \frac{1}{2} \xi \left(\xi - 1 \right) \\
     * N_2(r, s, t) &= 1 - \xi^2 \\
     * N_3(r, s, t) &= \frac{1}{2} \xi \left(\xi + 1 \right)
     * \f}
     */
    void precompute_shape_functions();
};
}
