
#pragma once

#include "LinearSolver.hpp"

// Mumps includes
#include <MUMPS/dmumps_c.h>
#include <MUMPS/smumps_c.h>

namespace neon
{
template <typename T>
struct MUMPSWrapper;

template <>
struct MUMPSWrapper<float>
{
    using MUMPS_STRUC_C = SMUMPS_STRUC_C;
    static void mumps_c(MUMPS_STRUC_C& info) { smumps_c(&info); }
};

template <>
struct MUMPSWrapper<double>
{
    using MUMPS_STRUC_C = DMUMPS_STRUC_C;
    static void mumps_c(MUMPS_STRUC_C& info) { dmumps_c(&info); }
};

/**
 * MUMPS is a base class for the multifrontal direct solver.  This solver
 * is widely used in the parallel solution of linear systems.
 *
 * TODO Put link to the solver website and documentation
 */
class MUMPS : public DirectLinearSolver
{
public:
    enum Ordering { AMD, AMF = 2, Scotch, Pord, Metis, QAMD, Automatic };

    // Jobs in MUMPS use the following:
    // 4   Job = 1 && Job = 2
    // 5   Job = 2 && Job = 3
    // 6   Job = 1 && Job = 2 && Job = 3
    enum Job {
        Terminate = -2,
        Initialisation = -1,
        Analysis = 1,
        Factorisation = 2,
        BackSubstitution = 3
    };

    enum Residual { None, Expensive, Cheap };

    enum MatrixProperty { Unsymmetric, SPD, GeneralSymmetric };

    using MUMPSAdapter = MUMPSWrapper<SparseMatrix::Scalar>;

public:
    explicit MUMPS(MatrixProperty const symmetric_flag);

    ~MUMPS();

protected:
    /**
     * Expand the sparse matrix into coordinate format only using the upper
     * diagonal values if the only_upper flag is set, otherwise expand
     * the entire matrix into compressed coordinate (COO) format
     */
    virtual void allocate_coordinate_format_storage(SparseMatrix const& A) = 0;

    void internal_solve(SparseMatrix const& A, vector& x, vector const& b);

protected:
    MUMPSAdapter::MUMPS_STRUC_C info;

    std::vector<int> rows, cols;      //!< Row and column index storage (uncompressed)
    std::vector<double> coefficients; //!< Sparse matrix coefficients
};

/**
 * MUMPSLLT is the LL^T factorisation (Cholesky) for a symmetric positive
 * definite matrix.  This solver can only be applied on a linear system and
 * takes the lower triangular part of the sparse matrix
 */
class MUMPSLLT : public MUMPS
{
public:
    MUMPSLLT() : MUMPS(MUMPS::MatrixProperty::SPD) {}

    void solve(SparseMatrix const& A, vector& x, vector const& b) override final;

protected:
    virtual void allocate_coordinate_format_storage(SparseMatrix const& A) override final;
};

/**
 * MUMPSLLT is the LL^T factorisation (Cholesky) for a general unsymmetric matrix.
 * This solver can only be applied on a linear system and takes the entire
 * uncompressed matrix
 */
class MUMPSLU : public MUMPS
{
public:
    MUMPSLU() : MUMPS(MUMPS::MatrixProperty::Unsymmetric) {}

    void solve(SparseMatrix const& A, vector& x, vector const& b) override final;

protected:
    virtual void allocate_coordinate_format_storage(SparseMatrix const& A) override final;
};
}
