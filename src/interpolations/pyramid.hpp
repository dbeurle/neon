
#pragma once

#include "shape_function.hpp"

namespace neon
{
/// pyramid5 is responsible for computing the shape functions for an prism or
/// wedge shaped six noded element.
/// The shape functions and nodal orderings are taken from \cite Bedrosian1992
/// with a volume of the unit element of one
/**
 * Initialise the shape functions to the following polynomials
 * \f{align*}{
 * N_1(r, s, t) &= \frac{1}{4}\left[(1+r)*(1+s) - t + \frac{rst}{1-t}\right] \\
 * N_2(r, s, t) &= \frac{1}{4}\left[(1-r)*(1+s) - t - \frac{rst}{1-t}\right] \\
 * N_3(r, s, t) &= \frac{1}{4}\left[(1-r)*(1-s) - t + \frac{rst}{1-t}\right] \\
 * N_4(r, s, t) &= \frac{1}{4}\left[(1+r)*(1-s) - t - \frac{rst}{1-t}\right] \\
 * N_5(r, s, t) &= t
 * \f}
 */
class pyramid5 : public volume_interpolation
{
public:
    explicit pyramid5();

    /// Evaluate the shape functions at the natural coordinate
    virtual auto evaluate(coordinate_type const& coordinate) const noexcept(false)
        -> value_type override final;
};

/// pyramid13 is responsible for computing the shape functions for a pyramid
/// twelve noded element. The shape functions and nodal orderings are taken from
/// \cite Bedrosian1992 with a volume of the unit element of one
class pyramid13 : public volume_interpolation
{
public:
    explicit pyramid13();

    /// Evaluate the shape functions at the natural coordinate
    virtual auto evaluate(coordinate_type const& coordinate) const noexcept(false)
        -> value_type override final;
};
}
