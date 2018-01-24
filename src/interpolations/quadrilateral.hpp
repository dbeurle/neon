
#pragma once

#include "shape_function.hpp"
#include "quadrature/quadrilateral_quadrature.hpp"

namespace neon
{
/** A finite element with 4 nodal points and an isoparametric formulation */
class quadrilateral4 : public surface_interpolation
{
public:
    quadrilateral4(quadrilateral_quadrature::Rule rule);

    virtual int nodes() const override final { return 4; }

    [[nodiscard]] double compute_measure(matrix const& nodal_coordinates) const;

protected:
    /**
     * Initialize the shape functions of the quadrilateral 4
     * node element to be:
     * \f{align*}{
     * N_1(\xi, \eta) &= \frac{1}{4}(1-\xi)(1-\eta) \\
     * N_2(\xi, \eta) &= \frac{1}{4}(1+\xi)(1-\eta) \\
     * N_3(\xi, \eta) &= \frac{1}{4}(1+\xi)(1+\eta) \\
     * N_4(\xi, \eta) &= \frac{1}{4}(1-\xi)(1+\eta)
     * \f}
     * TODO The derivatives of the surface element
     */
    void precompute_shape_functions();
};

/** A finite element with 8 nodal points and an isoparametric formulation */
class quadrilateral8 : public surface_interpolation
{
public:
    quadrilateral8(quadrilateral_quadrature::Rule rule);

    int nodes() const override final { return 8; }

    [[nodiscard]] double compute_measure(matrix const& nodal_coordinates) const;

protected:
    void precompute_shape_functions();
};

/** A finite element with 8 nodal points and an isoparametric formulation */
class quadrilateral9 : public surface_interpolation
{
public:
    quadrilateral9(quadrilateral_quadrature::Rule rule);

    int nodes() const override final { return 9; }

    [[nodiscard]] double compute_measure(matrix const& nodal_coordinates) const;

protected:
    void precompute_shape_functions();
};
}
