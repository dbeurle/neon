
#pragma once

#include "shape_function.hpp"
#include "quadrature/pyramid_quadrature.hpp"

namespace neon
{
/// pyramid5 is responsible for computing the shape functions for an prism or
/// wedge shaped six noded element.
/// The shape functions and nodal orderings are taken from \cite Bedrosian1992
/// with a volume of the unit element of one
class pyramid5 : public volume_interpolation
{
public:
    explicit pyramid5(pyramid_quadrature::point const p);

    virtual ~pyramid5() = default;

    int nodes() const override final { return 5; }

    double compute_measure(matrix3x const& nodal_coordinates) const;

protected:
    /**
     * Initialise the shape functions to the following polynomials
     * \f{align*}{
     * N_1(r, s, t) &= \frac{1}{4}\left[(1+r)*(1+s) - t + \frac{rst}{1-t}\right] \\
     * N_2(r, s, t) &= \frac{1}{4}\left[(1-r)*(1+s) - t - \frac{rst}{1-t}\right] \\
     * N_3(r, s, t) &= \frac{1}{4}\left[(1-r)*(1-s) - t + \frac{rst}{1-t}\right] \\
     * N_4(r, s, t) &= \frac{1}{4}\left[(1+r)*(1-s) - t - \frac{rst}{1-t}\right] \\
     * N_5(r, s, t) &= t \\
     * \f}
     */
    void precompute_shape_functions();
};

/// pyramid12 is responsible for computing the shape functions for a pyramid
/// twelve noded element. The shape functions and nodal orderings are taken from
/// \cite Bedrosian1992 with a volume of the unit element of one
class pyramid12 : public volume_interpolation
{
public:
    explicit pyramid12(pyramid_quadrature::point const p);

    virtual ~pyramid12() = default;

    int nodes() const override final { return 12; }

    double compute_measure(matrix3x const& nodal_coordinates) const;

protected:
    /**
     * Initialise the shape functions to the following:
     * \f{align*}{
     * N_1(r, s, t, \xi) &= r(2r - 1) \frac{1}{2} \xi (\xi - 1) \\
     * N_2(r, s, t, \xi) &= s(2s - 1) \frac{1}{2} \xi (\xi - 1) \\
     * N_3(r, s, t, \xi) &= t(2t - 1) \frac{1}{2} \xi (\xi - 1) \\
     * N_4(r, s, t, \xi) &= 4rs \frac{1}{2} \xi (\xi - 1) \\
     * N_5(r, s, t, \xi) &= 4st \frac{1}{2} \xi (\xi - 1) \\
     * N_6(r, s, t, \xi) &= 4rt \frac{1}{2} \xi (\xi - 1) \\
     * N_7(r, s, t, \xi) &= r (1 - \xi^2) \\
     * N_8(r, s, t, \xi) &= s (1 - \xi^2) \\
     * N_9(r, s, t, \xi) &= t (1 - \xi^2) \\
     * N_10(r, s, t, \xi) &= r(2r - 1) \frac{1}{2} \xi (\xi + 1) \\
     * N_11(r, s, t, \xi) &= s(2s - 1) \frac{1}{2} \xi (\xi + 1) \\
     * N_12(r, s, t, \xi) &= t(2t - 1) \frac{1}{2} \xi (\xi + 1) \\
     * N_13(r, s, t, \xi) &= 4rs \frac{1}{2} \xi (\xi + 1) \\
     * N_14(r, s, t, \xi) &= 4st \frac{1}{2} \xi (\xi + 1) \\
     * N_15(r, s, t, \xi) &= 4rt \frac{1}{2} \xi (\xi + 1) \\
     * \f}
     * where
     * \f{align*}
     *     t = 1 - r - s
     * \f}
     */
    void precompute_shape_functions();
};
}
