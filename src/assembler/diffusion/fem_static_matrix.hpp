
#pragma once

#include "io/file_output.hpp"
#include "mesh/diffusion/fem_mesh.hpp"
#include "numeric/sparse_matrix.hpp"

namespace neon
{
class linear_solver;
}

namespace neon::diffusion
{
/// fem_static_matrix is responsible for formulating and solving the steady state
/// scalar diffusion equation.
///   \f{align*}{
///    \mathbf{K} \mathbf{d} &= \mathbf{f}
///   \f}
/// This performs the gather operations for the boundary
/// conditions and applies the Dirichlet conditions to the matrix.  For dynamic
/// problems see fem_dynamic_matrix.
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

    /// Cache flag for matrix assembly
    bool is_sparsity_computed{false};
    /// Conductivity matrix
    sparse_matrix K;
    /// Heat vector
    vector f;
    /// Temperature vector
    vector d;

    std::unique_ptr<linear_solver> solver;
};
}
