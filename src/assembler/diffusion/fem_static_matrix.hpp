
#pragma once

#include "io/FileIO.hpp"
#include "mesh/diffusion/fem_mesh.hpp"
#include "numeric/sparse_matrix.hpp"

namespace neon
{
class LinearSolver;
}

namespace neon::diffusion
{
/**
 * fem_static_matrix is responsible for formulating and solving the steady state
 * scalar diffusion equation.
   \f{align*}{
    \mathbf{K} \mathbf{d} &= \mathbf{f}
   \f}
 * This performs the gather operations for the boundary
 * conditions and applies the Dirichlet conditions to the matrix.  For dynamic
 * problems see fem_dynamic_matrix.
 *
 */
class fem_static_matrix
{
public:
    explicit fem_static_matrix(fem_mesh& mesh, json const& simulation_data);

    ~fem_static_matrix();

    /** Computes the steady state solution of the diffusion equation */
    virtual void solve();

protected:
    /**
     * Compute the sparse pattern of the coefficient matrix but use the doublet
     * list as a means of forming the sparsity pattern.  This is a memory
     * intensive operation and should be replaced in the future
     */
    void compute_sparsity_pattern();

    /** Assembles the external contribution vector */
    void compute_external_force(double const load_factor = 1.0);

    /** Assembles the conductivity matrix */
    void assemble_stiffness();

    /**
     * Apply dirichlet conditions to the system defined by A, x, and b.
     * This method selects the row and column of the degree of freedom with an
     * imposed Dirichlet condition.
     *
     * The algorithm takes the row and performs a zeroing operation on it.  To
     * satisfy the equation system, the constrained column is multiplied with
     * the Dirichlet value and subtracted from the matching DoF on the right
     * hand side vector.  The column is then zeroed and the diagonal DoF is
     * then corrected such that \f$ A_{dof} * x_{dof} = f_{dof} == A_{dof} * x_{dof} \f$ so the
     * equation system is satisfied.
     *
     * For inner and outer vector reference see
     * https://eigen.tuxfamily.org/dox/group__TutorialSparse.html
     */
    void apply_dirichlet_conditions(sparse_matrix& A, vector& x, vector& b);

protected:
    fem_mesh& mesh;

    bool is_sparsity_computed = false;

    sparse_matrix K; /// Conductivity matrix
    vector f;        /// Heat vector
    vector d;        /// Temperature vector

    FileIO file_io;

    std::unique_ptr<LinearSolver> linear_solver;
};
}
