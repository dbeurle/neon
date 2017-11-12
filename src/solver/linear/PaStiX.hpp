
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

    bool sparsity_pattern_changed = true;
};

/**
 * PaStiXLU is a supernodal direct solver with multithreading support.  This
 * performs the LU factorisation of a square matrix system
 */
class PaStiXLU : public DirectLinearSolver
{
public:
    void solve(SparseMatrix const& A, Vector& x, Vector const& b) override final;
};
}
