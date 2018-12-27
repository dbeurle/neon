
#pragma once

#include "shape_function.hpp"

namespace neon
{
/// Triangular 3 node element (constant-strain / stress element).  Usually
/// recommended for solving the Laplace equation but exhibits poor convergence
/// and is generally not recommended for production use.
class triangle3 : public surface_interpolation
{
public:
    triangle3();

    /// Evaluate the shape functions at the natural coordinate
    virtual auto evaluate(coordinate_type const& coordinate) const noexcept(false)
        -> value_type override final;
};

/// A six node triangle with midside nodes and quadrature shape functions.
/// \sa Hughes2012
class triangle6 : public surface_interpolation
{
public:
    triangle6();

    /// Evaluate the shape functions at the natural coordinate
    virtual auto evaluate(coordinate_type const& coordinate) const noexcept(false)
        -> value_type override final;
};
}
