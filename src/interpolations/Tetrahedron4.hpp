
#pragma once

#include "ShapeFunction.hpp"

namespace neon
{
/**
 * Tetrahedron4 is an isoparametric tetrahedral 4 node element with analytical
 * integration.
 */
class Tetrahedron4 : public VolumeInterpolation
{
public:
    Tetrahedron4();

protected:
    /**
     * Compute six times the volume of the tetrahedral element
     * @param nodal coordinates
     * @return 6 * volume of the element
     */
    double computeSixV(Matrix const& nodalCoordinates) const;

    /**
     * Initialize the derivative matrix.  Rows are the derivatives w.r.t the
     * tetrahedral coordinates and the columns are the shape functions (N_1 to N_4)
     */
    void setDerivativeMatrix();
    void initializeSource();
    void initializeMass();

    Matrix dN, mass;
    Vector force;
};
}
