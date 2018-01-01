
#pragma once

#include "LinearSolver.hpp"

#include <Eigen/Sparse>

#include <Eigen/PaStiXSupport>

// Respect this order for access functions from libpastix because of silly design
extern "C" {
#include "pastix.h"

#include "cscd_utils.h"

#include "read_matrix.h"

#include "get_options.h"
}

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

    void solve(SparseMatrix const& A, vector& x, vector const& b) override final;

private:
    Eigen::PastixLDLT<Eigen::SparseMatrix<double>, Eigen::Upper> ldlt;

    std::array<pastix_int_t, IPARM_SIZE> iparm; /* integer parameters for pastix */
    std::array<double, DPARM_SIZE> dparm;       /* floating parameters for pastix */
};

/**
 * PaStiXLU is a supernodal direct solver with multithreading support.  This
 * performs the LU factorisation of a square matrix system
 */
class PaStiXLU : public DirectLinearSolver
{
public:
    PaStiXLU();

    void solve(SparseMatrix const& A, vector& x, vector const& b) override final;

private:
    // BUG Likely not going to work with unsymmetric matrix because of row and
    // column ordering change.  Should give the transpose of the matrix but
    // unsure why Eigen can't handle this
    Eigen::PastixLU<Eigen::SparseMatrix<double>> lu;
};
}
