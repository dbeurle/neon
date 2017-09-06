
#include "MUMPS.hpp"
#include "Exceptions.hpp"

#include <Eigen/Sparse>

// Mumps includes
#include <MUMPS/dmumps_c.h>
#include <MUMPS/smumps_c.h>

#include <chrono>
#include <iostream>

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

void MUMPSLU::allocate_coordinate_format_storage(SparseMatrix const& A)
{
    rows.clear();
    cols.clear();
    coefficients.clear();

    rows.reserve(A.nonZeros());
    cols.reserve(A.nonZeros());
    coefficients.reserve(A.nonZeros());

    for (auto k = 0, nonzeros = 0; k < A.outerSize(); ++k)
    {
        for (SparseMatrix::InnerIterator it(A, k); it; ++it)
        {
            coefficients.emplace_back(it.value());
            rows.emplace_back(it.row() + 1);
            cols.emplace_back(it.col() + 1);
        }
    }
}

void MUMPSLLT::allocate_coordinate_format_storage(SparseMatrix const& A)
{
    rows.clear();
    cols.clear();
    coefficients.clear();

    rows.reserve(A.nonZeros());
    cols.reserve(A.nonZeros());
    coefficients.reserve(A.nonZeros());

    // Decompress the upper part of the sparse matrix
    for (auto k = 0; k < A.outerSize(); ++k)
    {
        for (SparseMatrix::InnerIterator it(A, k); it; ++it)
        {
            if (it.col() >= it.row())
            {
                coefficients.emplace_back(it.value());
                rows.emplace_back(it.row() + 1);
                cols.emplace_back(it.col() + 1);
            }
        }
    }
}

void MUMPS::internal_solve(SparseMatrix const& A, Vector& x, Vector const& b, int const symmetric_flag)
{
    auto start = std::chrono::high_resolution_clock::now();

    x = b;

    using MUMPSAdapter = MUMPSWrapper<double>;

    MUMPSAdapter::MUMPS_STRUC_C info;

    info.job = Job::Initialisation;

    // Par determines parallel data location
    // 0  : Host is not involved with computes
    // 1  : Host is involved
    info.par = 1;
    info.comm_fortran = -987654;

    // Symmetric matrix flag
    info.sym = symmetric_flag;

    MUMPSAdapter::mumps_c(info);

    // Verbosity levels
    info.icntl[0] = 1;
    info.icntl[1] = 1;
    info.icntl[2] = 1;
    info.icntl[3] = 1;

    // Ordering algorithm
    info.icntl[6] = Ordering::Automatic;

    // Iterative refinement
    info.icntl[9] = 100;

    // Compute the residual
    info.icntl[10] = Residual::Cheap;

    info.n = A.rows();
    info.nz = coefficients.size();

    info.a = coefficients.data();
    info.irn = rows.data();
    info.jcn = cols.data();

    // Analysis phase
    info.job = 1;
    MUMPSAdapter::mumps_c(info);

    if (info.info[0] < 0) throw computational_error("Error in analysis phase of MUMPS solver\n");

    // Factorization phase
    info.job = Job::Factorisation;
    MUMPSAdapter::mumps_c(info);

    if (info.info[0] < 0)
        throw computational_error("Error in factorisation phase of MUMPS solver\n");

    info.rhs = x.data();
    info.nrhs = 1;
    info.lrhs = info.n;

    info.job = Job::BackSubstitution;
    MUMPSAdapter::mumps_c(info);

    if (info.info[0] < 0)
        throw computational_error("Error in back substitution phase of MUMPS solver\n");

    // Take out the trash
    info.job = Job::Terminate;
    MUMPSAdapter::mumps_c(info);

    if (info.info[0] < 0) throw computational_error("Error in cleanup phase\n");

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;

    std::cout << std::string(6, ' ') << "MUMPS solver took " << elapsed_seconds.count() << "s\n";
}

void MUMPSLU::solve(SparseMatrix const& A, Vector& x, Vector const& b)
{
    this->allocate_coordinate_format_storage(A);
    internal_solve(A, x, b, MatrixProperty::Unsymmetric);
}

void MUMPSLLT::solve(SparseMatrix const& A, Vector& x, Vector const& b)
{
    this->allocate_coordinate_format_storage(A);
    internal_solve(A, x, b, MatrixProperty::SPD);
}
}
