
#pragma once

#include "io/FileIO.hpp"
#include "mesh/diffusion/femMesh.hpp"
#include "numeric/SparseTypes.hpp"

#include <json/forwards.h>

namespace neon
{
class LinearSolver;
}

namespace neon::diffusion
{
/**
 * femStaticMatrix is responsible for formulating and solving the steady state
 * scalar diffusion equation.
   \f{align*}{
    \mathbf{K} \mathbf{d} &= \mathbf{f}
   \f}
 * This performs the gather operations for the boundary
 * conditions and applies the Dirichlet conditions to the matrix.  For dynamic
 * problems see femDynamicMatrix.
 *
 *
 */
class femStaticMatrix
{
public:
    explicit femStaticMatrix(femMesh& fem_mesh, Json::Value const& solver_data, FileIO&& file_io);

    ~femStaticMatrix();

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
    void apply_dirichlet_conditions(SparseMatrix& A, Vector& x, Vector& b);

protected:
    femMesh& fem_mesh;

    bool is_sparsity_computed = false;

    SparseMatrix K; //!< Conductivity matrix
    Vector f;       //!< Heat vector
    Vector d;       //!< Temperature vector

    FileIO file_io;

    std::unique_ptr<LinearSolver> linear_solver;
};
}
