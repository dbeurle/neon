
#pragma once

#include "ShapeFunction.hpp"
#include "quadrature/HexahedronQuadrature.hpp"

namespace neon
{
class Hexahedron8 : public VolumeInterpolation
{
public:
    Hexahedron8(HexahedronQuadrature::Rule rule);

    virtual int nodes() const override final { return 8; }

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
}
