
#pragma once

/// @file

#include "numeric/dense_matrix.hpp"
#include "numeric/sparse_matrix.hpp"

namespace neon
{
/// linear_solver is to setup a linear solver with designated parameters from
/// the input file.  This is the interface for every linear solver in neon
class linear_solver
{
public:
    virtual ~linear_solver() = default;

    virtual void solve(sparse_matrix const& A, vector& x, vector const& b) = 0;

    /// Notifies the linear solvers of a change in sparsity structure of A
    void update_sparsity_pattern() { build_sparsity_pattern = true; }

protected:
    bool build_sparsity_pattern{true};
};

class iterative_linear_solver : public linear_solver
{
public:
    /// Construct with default for residual and maximum iterations
    explicit iterative_linear_solver() = default;

    /// Override the default for the residual tolerance
    explicit iterative_linear_solver(double const residual_tolerance);

    /// Override the default for the maximum number of iterations
    explicit iterative_linear_solver(std::int32_t const max_iterations);

    /// Override the default for the residual tolerance and maximum number of iterations
    explicit iterative_linear_solver(double const residual_tolerance,
                                     std::int32_t const max_iterations);

protected:
    void compute_symmetric_reordering(sparse_matrix const& input_matrix);

    void apply_permutation(sparse_matrix const& input_matrix, vector const& input_rhs);

protected:
    double residual_tolerance{1.0e-5};
    std::int32_t max_iterations{2000};

    sparse_matrix A;
    vector b;

    permutation_matrix P;
};

/// conjugate_gradient is a simple solver wrapper for the preconditioned conjugate gradient
/// solver from Eigen.  This is a multithreaded solver when beneficial.
/// The preconditioner available is incomplete Cholesky, LU and simple Jacobi.
///
/// The benefit of this solver is the ability to use a previous
/// solution as a starting point.  This is useful in time analyses
/// when the solution is not expected to change significantly.
class conjugate_gradient : public iterative_linear_solver
{
public:
    using iterative_linear_solver::iterative_linear_solver;

    void solve(sparse_matrix const& input_matrix, vector& x, vector const& input_rhs) override final;
};

/// biconjugate_gradient_stabilised is a simple solver wrapper for the preconditioned bi-conjugate gradient
/// stabilised solver from Eigen.  This is multithreaded when beneficial.
/// The preconditioner available is incomplete Cholesky, LU and simple Jacobi.
///
/// The benefit of this solver is the ability to use a previous
/// solution as a starting point.  This is useful in time analyses
/// as the solution is not expected to change significantly.
class biconjugate_gradient_stabilised : public iterative_linear_solver
{
public:
    using iterative_linear_solver::iterative_linear_solver;

    void solve(sparse_matrix const& A, vector& x, vector const& b) override final;
};

class direct_linear_solver : public linear_solver
{
};

/// SparseLU is a single threaded sparse LU factorization using AMD reordering.
/// This solver is not recommended over the industrial grade solver PaStiX when
/// using a direct solver except for small problems or when PaStiX is not available
class SparseLU : public direct_linear_solver
{
public:
    void solve(sparse_matrix const& A, vector& x, vector const& b) override final;

private:
    Eigen::SparseLU<sparse_matrix, Eigen::AMDOrdering<std::int32_t>> lu;
};

/// SparseLLT is a single threaded sparse Cholesky factorization using AMD reordering.
/// This solver is not recommended over the industrial grade solver PaStiX when
/// using a direct solver except for small problems or when PaStiX is not available
class SparseLLT : public direct_linear_solver
{
public:
    void solve(sparse_matrix const& A, vector& x, vector const& b) override final;

private:
    Eigen::SimplicialLLT<Eigen::SparseMatrix<sparse_matrix::Scalar>> llt;
};
}
