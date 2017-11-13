
#pragma once

#include "ShapeFunction.hpp"
#include "quadrature/HexahedronQuadrature.hpp"

namespace neon
{
/**
 * Hexahedron8 is responsible for computing the tri-linear shape functions for an
 * eight noded hexahedron element.
 * The shape functions and ordering is taken from Hughes 1986 - Linear static and dynamic finite
 * elements.
 */
class Hexahedron8 : public VolumeInterpolation
{
public:
    explicit Hexahedron8(HexahedronQuadrature::Rule rule);

    int nodes() const override final { return 8; }

    double compute_measure(Matrix const& nodal_coordinates) const;

protected:
    /**
     * Initialize the shape functions to the following polynomials
     * \f{align*}{
     * N_1(\xi, \eta, \zeta) &= \frac{1}{8}(1-\xi)(1-\eta)(1-\zeta) \\
     * N_2(\xi, \eta, \zeta) &= \frac{1}{8}(1+\xi)(1-\eta)(1-\zeta) \\
     * N_3(\xi, \eta, \zeta) &= \frac{1}{8}(1+\xi)(1+\eta)(1-\zeta) \\
     * N_4(\xi, \eta, \zeta) &= \frac{1}{8}(1-\xi)(1+\eta)(1-\zeta) \\
     * N_5(\xi, \eta, \zeta) &= \frac{1}{8}(1-\xi)(1-\eta)(1+\zeta) \\
     * N_6(\xi, \eta, \zeta) &= \frac{1}{8}(1+\xi)(1-\eta)(1+\zeta) \\
     * N_7(\xi, \eta, \zeta) &= \frac{1}{8}(1+\xi)(1+\eta)(1+\zeta) \\
     * N_8(\xi, \eta, \zeta) &= \frac{1}{8}(1-\xi)(1+\eta)(1+\zeta)
     * \f}
     */
    void precompute_shape_functions();
};

/**
 * Hexahedron20 is responsible for computing the quadratic shape functions for an
 * twenty noded hexahedron element.  Nodes are only defined on the midside
 * and corner nodes.  The node ordering is from Hughes.
 */
class Hexahedron20 : public VolumeInterpolation
{
public:
    explicit Hexahedron20(HexahedronQuadrature::Rule rule);

    int nodes() const override final { return 20; }

    double compute_measure(Matrix const& nodal_coordinates) const;

protected:
    void precompute_shape_functions();
};

/**
 * Hexahedron27 is responsible for computing the quadratic shape functions for an
 * twenty-nine noded hexahedron element.  Nodes are also on the faces and the centre
 * of the reference cube.  The node ordering is from Hughes.
 */
class Hexahedron27 : public VolumeInterpolation
{
public:
    explicit Hexahedron27(HexahedronQuadrature::Rule rule);

    int nodes() const override final { return 27; }

    double compute_measure(Matrix const& nodal_coordinates) const;

protected:
    void precompute_shape_functions();
};
}
