
#pragma once

#include "numeric/DenseTypes.hpp"
#include "numeric/SparseTypes.hpp"

namespace neon
{
/**
 * LinearSolver is to setup a linear solver with designated parameters from
 * the input file.  This is the interface for every linear solver in neon
 */
class LinearSolver
{
public:
    virtual void solve(SparseMatrix const& A, Vector& x, Vector const& b) = 0;

    /** Notifies the linear solvers to of a change in sparsity structure of A */
    void update_sparsity_pattern() { build_sparsity_pattern = true; }

protected:
    bool build_sparsity_pattern = true;
};

class IterativeLinearSolver : public LinearSolver
{
public:
    explicit IterativeLinearSolver() = default;
    explicit IterativeLinearSolver(double const residual_tolerance);
    explicit IterativeLinearSolver(int const max_iterations);
    explicit IterativeLinearSolver(double const residual_tolerance, int const max_iterations);

protected:
    double residual_tolerance = 1.0e-5;
    int max_iterations = 1000;
};

/**
 * ConjugateGradient is a simple solver wrapper for the preconditioned conjugate gradient
 * solver from Eigen.  This is a multithreaded solver when beneficial.
 * The preconditioner available is incomplete Cholesky, LU and simple Jacobi.
 *
 * The benefit of this solver is the ability to use a previous
 * solution as a starting point.  This is useful in time analyses
 * when the solution is not expected to change significantly.
 */
class ConjugateGradient : public IterativeLinearSolver
{
public:
    using IterativeLinearSolver::IterativeLinearSolver;

    void solve(SparseMatrix const& A, Vector& x, Vector const& b) override final;
};

/**
 * BiCGStab is a simple solver wrapper for the preconditioned bi-conjugate gradient
 * stabilised solver from Eigen.  This is a multithreaded solver when beneficial.
 * The preconditioner available is incomplete Cholesky, LU and simple Jacobi.
 *
 * The benefit of this solver is the ability to use a previous
 * solution as a starting point.  This is useful in time analyses
 * as the solution is not expected to change significantly.
 */
class BiCGStab : public IterativeLinearSolver
{
public:
    using IterativeLinearSolver::IterativeLinearSolver;

    void solve(SparseMatrix const& A, Vector& x, Vector const& b) override final;
};

class DirectLinearSolver : public LinearSolver
{
};

/**
 * SparseLU is a single threaded sparse LU factorization using AMD reordering.
 * This solver is not recommended over the industrial grade solver PaStiX when
 * using a direct solver except for small problems or when PaStiX is not available
 */
class SparseLU : public DirectLinearSolver
{
public:
    void solve(SparseMatrix const& A, Vector& x, Vector const& b) override final;
};

/**
 * SparseLLT is a single threaded sparse Cholesky factorization using AMD reordering.
 * This solver is not recommended over the industrial grade solver PaStiX when
 * using a direct solver except for small problems or when PaStiX is not available
 */
class SparseLLT : public DirectLinearSolver
{
public:
    void solve(SparseMatrix const& A, Vector& x, Vector const& b) override final;
};
}
