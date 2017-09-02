
#pragma once

#include "LinearSolver.hpp"

namespace neon
{
/**
 * PaStiXLDLT is a supernodal direct solver with multithreading support.  This
 * performs the LDLT Cholesky factorisation of a symmetric system
 */
class PaStiXLDLT : public DirectLinearSolver
{
public:
    void solve(SparseMatrix const& A, Vector& x, Vector const& b) override final;
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
