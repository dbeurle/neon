
#pragma once

#include "solver/linear/linear_solver.hpp"

namespace neon::fem::detail
{
/// fem_linear_static_matrix is responsible for assembling the stiffness matrix,
/// solving the linear system of equations and dispatching the postprocessor
/// for a mesh with linear (or infinitesimal) behaviour.  This applies to
/// any mesh that provides the following methods:
/// \p tangent_stiffness
/// \p dirichlet_boundaries
/// and these are used to form the linear system \p K*d=f_ext and this is then
/// passed through the linear solver and the postprocessor.
template <typename fem_mesh_type>
class linear_static_matrix
{
public:
    using mesh_type = fem_mesh_type;

public:
    linear_static_matrix(mesh_type& mesh) : mesh(mesh) {}

    ~linear_static_matrix() = default;

    void solve() {}

protected:
    /// Compute the sparse pattern of the coefficient matrix but use the doublet
    /// list as a means of forming the sparsity pattern.  This is a memory
    /// intensive operation and should be replaced in the future
    void compute_sparsity_pattern() {}

    /// Assembles the external contribution vector
    void compute_external_force(double const load_factor = 1.0) {}

    /// Assembles the linear stiffness matrix
    void assemble_stiffness() {}

protected:
    mesh_type& mesh;

    std::unique_ptr<LinearSolver> linear_solver;

    /// Stiffness matrix
    sparse_matrix K;

    /// Right hand side load vector
    vector f_ext;

    /// Solution vector
    vector d;
};
}
