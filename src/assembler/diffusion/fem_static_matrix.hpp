
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
 */
class fem_static_matrix
{
public:
    explicit fem_static_matrix(fem_mesh& mesh, json const& simulation_data);

    ~fem_static_matrix();

    /// Computes the steady state solution of the diffusion equation
    virtual void solve();

protected:
    /// Compute the sparse pattern of the coefficient matrix but use the doublet
    /// list as a means of forming the sparsity pattern.  This is a memory
    /// intensive operation and should be replaced in the future
    void compute_sparsity_pattern();

    /// Assembles the external contribution vector
    void compute_external_force(double const load_factor = 1.0);

    /// Assembles the conductivity matrix
    void assemble_stiffness();

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
