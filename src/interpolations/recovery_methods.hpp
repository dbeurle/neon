
#pragma once

#include "numeric/dense_matrix.hpp"

namespace neon
{
class recovery_method
{
};

class local_extrapolation : public recovery_method
{
public:
    /// Compute the extrapolation matrix to allow for quadrature valued variables
    /// to be averaged to the nodal points without ill-effects when using a
    /// least squares (for example with quadratric tetrahedron elements)
    /// developed in \cite Durand2014
    void compute_extrapolation_matrix(matrix const& N,
                                      matrix const& local_nodal_coordinates,
                                      matrix const& local_quadrature_coordinates);

private:
    /// Quadrature point to nodal point mapping
    matrix extrapolation;
};

class super_convergent_patch_recovery : public recovery_method
{
};

}
