
#pragma once

/// @file

#include "linear_solver.hpp"

extern "C" {
#include <pastix.h>
#include <spm.h>
}

#include <memory>

namespace neon
{
class PaStiX : public direct_linear_solver
{
public:
    explicit PaStiX();

    virtual ~PaStiX();

protected:
    void analyse_pattern();

protected:
    std::array<double, DPARM_SIZE> m_float_parameters;
    std::array<pastix_int_t, IPARM_SIZE> m_integer_parameters;

    std::unique_ptr<spmatrix_t> m_matrix;

    pastix_data_t* m_data{nullptr};
};

/// PaStiXLDLT is a supernodal direct solver with multithreading support.  This
/// performs the LDLT Cholesky factorisation of a symmetric system
class PaStiXLDLT : public PaStiX
{
public:
    explicit PaStiXLDLT();

    void solve(sparse_matrix const& A, vector& x, vector const& b) override final;
};

/// PaStiXLU is a supernodal direct solver with multithreading support.  This
/// performs the LU factorisation of a square matrix system
class PaStiXLU : public PaStiX
{
public:
    explicit PaStiXLU();

    void solve(sparse_matrix const& A, vector& x, vector const& b) override final;
};
}
