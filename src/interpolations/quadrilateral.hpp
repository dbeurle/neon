
#pragma once

/// @file

#include "shape_function.hpp"

namespace neon
{
/// A finite element with 4 nodal points and an isoparametric formulation
/**
 * Initialize the shape functions of the quadrilateral 4
 * node element to be:
 * \f{align*}{
 * N_1(\xi, \eta) &= \frac{1}{4}(1-\xi)(1-\eta) \\
 * N_2(\xi, \eta) &= \frac{1}{4}(1+\xi)(1-\eta) \\
 * N_3(\xi, \eta) &= \frac{1}{4}(1+\xi)(1+\eta) \\
 * N_4(\xi, \eta) &= \frac{1}{4}(1-\xi)(1+\eta)
 * \f}
 */
class quadrilateral4 : public surface_interpolation
{
public:
    explicit quadrilateral4();

    /// Evaluate the shape functions at the natural coordinate
    virtual auto evaluate(coordinate_type const& coordinate) const noexcept(false)
        -> value_type override final;
};

/// A finite element with 8 nodal points and an isoparametric formulation
class quadrilateral8 : public surface_interpolation
{
public:
    explicit quadrilateral8();

    /// Evaluate the shape functions at the natural coordinate
    virtual auto evaluate(coordinate_type const& coordinate) const noexcept(false)
        -> value_type override final;
};

/// A finite element with 8 nodal points and an isoparametric formulation
class quadrilateral9 : public surface_interpolation
{
public:
    explicit quadrilateral9();

    /// Evaluate the shape functions at the natural coordinate
    virtual auto evaluate(coordinate_type const& coordinate) const noexcept(false)
        -> value_type override final;
};
}
