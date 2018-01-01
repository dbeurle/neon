
#pragma once

#include "distributed_linear_solver.hpp"

extern "C" {
#include "pastix.h"
}

namespace neon
{
/**
 * distributed_pastix_ldlt is a supernodal direct solver with multithreading
 * and domain decomposition.  This performs the LDLT Cholesky factorisation of
 * a symmetric system.
 */
class distributed_pastix_ldlt : public distributed_direct_linear_solver
{
public:
    distributed_pastix_ldlt();

    void solve(SparseMatrix const& A, vector& x, vector const& b) override final;

protected:
    void analyse_pattern();

private:
    std::array<pastix_int_t, IPARM_SIZE> integer_parameters; //!< integer parameters
    std::array<double, DPARM_SIZE> float_parameters;         //!< floating parameters

    pastix_data_t* pastix_data = nullptr; //!< Storage structure needed by pastix

    // Reverse permutation tabular
    std::vector<pastix_int_t> invp, perm;
};

/**
 * distributed_PaStiX_lu is a supernodal direct solver with multithreading
 * support.  This performs the LU factorisation of a square matrix system
 */
// class distributed_PaStiX_lu : public distributed_direct_linear_solver
// {
// public:
//     PaStiXLU();
//
//     void solve(SparseMatrix const& A, vector& x, vector const& b) override final;
//
// private:
//     // BUG Likely not going to work with unsymmetric matrix because of row and
//     // column ordering change.  Should give the transpose of the matrix but
//     // unsure why Eigen can't handle this
//     Eigen::PastixLU<Eigen::SparseMatrix<double>> lu;
// };
}
