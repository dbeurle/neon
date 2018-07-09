
#pragma once

#include "shape_function.hpp"
#include "quadrature/prism_quadrature.hpp"

namespace neon
{
/// prism6 is responsible for computing the shape functions for an prism or
/// wedge shaped six noded element.
/// The shape functions and nodal orderings are taken from \cite Hughes2012
/// with a volume of the unit element of one
class prism6 : public volume_interpolation
{
public:
    explicit prism6(prism_quadrature::point const p);

    virtual ~prism6() = default;

    int nodes() const override final { return 6; }

    double compute_measure(matrix3x const& nodal_coordinates) const;

protected:
    /**
     * Initialise the shape functions to the following polynomials
     * \f{align*}{
     * N_1(r, s, t, \xi) &= \frac{1}{2}(1-\xi)r \\
     * N_2(r, s, t, \xi) &= \frac{1}{2}(1-\xi)s \\
     * N_3(r, s, t, \xi) &= \frac{1}{2}(1-\xi)t \\
     * N_4(r, s, t, \xi) &= \frac{1}{2}(1+\xi)r \\
     * N_5(r, s, t, \xi) &= \frac{1}{2}(1+\xi)s \\
     * N_6(r, s, t, \xi) &= \frac{1}{2}(1+\xi)r \\
     * \f}
     */
    void precompute_shape_functions();
};

/// prism15 is responsible for computing the shape functions for an prism or
/// wedge shaped fifteen noded element.
/// The shape functions and nodal orderings are taken from \cite Hughes2012
/// with a volume of the unit element of one
class prism15 : public volume_interpolation
{
public:
    explicit prism15(prism_quadrature::point const p);

    virtual ~prism15() = default;

    int nodes() const override final { return 15; }

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
