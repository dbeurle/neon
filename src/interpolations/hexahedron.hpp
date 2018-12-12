
#pragma once

#include "shape_function.hpp"
#include "quadrature/hexahedron/hexahedron_quadrature.hpp"

namespace neon
{
/// hexahedron8 is responsible for computing the tri-linear shape functions for an
/// eight noded hexahedron element.
/// The shape functions and ordering is from \cite Hughes2012
class hexahedron8 : public volume_interpolation
{
public:
    explicit hexahedron8(hexahedron_quadrature::point const p);

    virtual ~hexahedron8() = default;

    int nodes() const override final { return 8; }

    /// \return The element volume
    double compute_measure(matrix const& nodal_coordinates) const;

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

/// hexahedron20 is responsible for computing the quadratic shape functions for an
/// twenty noded hexahedron element.  Nodes are only defined on the midside
/// and corner nodes.  The node ordering is from \cite Hughes2012.
class hexahedron20 : public volume_interpolation
{
public:
    explicit hexahedron20(hexahedron_quadrature::point const p);

    virtual ~hexahedron20() = default;

    int nodes() const override final { return 20; }

    /// \return The element volume
    double compute_measure(matrix const& nodal_coordinates) const;

protected:
    void precompute_shape_functions();
};

/// hexahedron27 is responsible for computing the quadratic shape functions for an
/// twenty-nine noded hexahedron element.  Nodes are also on the faces and the centre
/// of the reference cube.  The node ordering is from \cite Hughes2012.
class hexahedron27 : public volume_interpolation
{
public:
    explicit hexahedron27(hexahedron_quadrature::point const p);

    virtual ~hexahedron27() = default;

    int nodes() const override final { return 27; }

    /// \return The element volume
    double compute_measure(matrix const& nodal_coordinates) const;

protected:
    void precompute_shape_functions();
};
}
