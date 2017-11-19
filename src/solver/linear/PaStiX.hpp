
#pragma once

#include "LinearSolver.hpp"

#include <Eigen/Sparse>

#include <Eigen/PaStiXSupport>

namespace neon
{
/**
 * PaStiXLDLT is a supernodal direct solver with multithreading support.  This
 * performs the LDLT Cholesky factorisation of a symmetric system
 */
class PaStiXLDLT : public DirectLinearSolver
{
public:
    PaStiXLDLT();

    void solve(SparseMatrix const& A, Vector& x, Vector const& b) override final;

private:
    Eigen::PastixLDLT<Eigen::SparseMatrix<double>, Eigen::Upper> ldlt;
};

/**
 * PaStiXLU is a supernodal direct solver with multithreading support.  This
 * performs the LU factorisation of a square matrix system
 */
class PaStiXLU : public DirectLinearSolver
{
public:
    PaStiXLU();

    void solve(SparseMatrix const& A, Vector& x, Vector const& b) override final;

private:
    // BUG Likely not going to work with unsymmetric matrix because of row and
    // column ordering change.  Should give the transpose of the matrix but
    // unsure why Eigen can't handle this
    Eigen::PastixLU<Eigen::SparseMatrix<double>> lu;
};
}
