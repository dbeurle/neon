
#pragma once

#include "numeric/SparseTypes.hpp"
#include <string>

namespace neon
{
/**
 * LinearSolver is to setup a linear solver with designated parameters from
 * the input file.  This is the interface for every linear solver in neon
 */
class LinearSolver
{
public:
    LinearSolver();

    virtual void solve(SparseMatrix const& A, Vector& x, Vector const& b) = 0;

protected:
    /**
     * SolverParam is a lightweight struct to hold the name of the solver
     * and the desired tolerances and maximum iterations for iterative solvers.
     */
    struct SolverParam
    {
    public:
        SolverParam(double tol, int iter) : tolerance(tol), max_iterations(iter) {}
        double tolerance;
        int max_iterations;
        std::string name;
    };

protected:
    SolverParam solverParam; //!< POD for the solver parameters
};

/**
 * PaStiX is a supernodal direct solver with multithreading support
 */
class PaStiX : public LinearSolver
{
public:
    void solve(SparseMatrix const& A, Vector& x, Vector const& b) override final;
};

/**
 * MUMPS is a multifrontal direct solver
 */
class MUMPS : public LinearSolver
{
public:
    void solve(SparseMatrix const& A, Vector& x, Vector const& b) override final;
};

/**
 * SparseLU is a single threaded sparse LU factorization using AMD reordering.
 * This solver is not recommended over the industrial grade solver PaStiX when
 * using a direct solver except for small problems or when PaStiX is not available
 */
class SparseLU : public LinearSolver
{
public:
    void solve(SparseMatrix const& A, Vector& x, Vector const& b) override final;
};

/**
 * pCG is a simple solver wrapper for the preconditioned conjugate gradient
 * solver from Eigen.  This is a multithreaded solver when beneficial.
 * The preconditioner available is incomplete Cholesky, LU and simple Jacobi.
 *
 * The benefit of this solver is the ability to use a previous
 * solution as a starting point.  This is useful in time analyses
 * when the solution is not expected to change significantly.
 */
class pCG : public LinearSolver
{
public:
    pCG() = default;
    pCG(double tol);
    pCG(int maxIter);
    pCG(double tol, int maxIter);

    void solve(SparseMatrix const& A, Vector& x, Vector const& b) override final;
};

/**
 * BiCGSTAB is a simple solver wrapper for the preconditioned bi-conjugate gradient
 * stabilized solver from Eigen.  This is a multithreaded solver when beneficial.
 * The preconditioner available is incomplete Cholesky, LU and simple Jacobi.
 *
 * The benefit of this solver is the ability to use a previous
 * solution as a starting point.  This is useful in time analyses
 * as the solution is not expected to change significantly.
 */
class BiCGSTAB : public LinearSolver
{
public:
    BiCGSTAB() = default;
    BiCGSTAB(double tol);
    BiCGSTAB(int maxIter);
    BiCGSTAB(double tol, int maxIter);

    void solve(SparseMatrix const& A, Vector& x, Vector const& b) override final;
};
}
