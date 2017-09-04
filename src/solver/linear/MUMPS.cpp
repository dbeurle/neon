
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

void MUMPS::allocate_coordinate_format_storage(SparseMatrix const& A, const bool only_upper)
{
    // Temporary matrix storage
    if (only_upper)
    {
        irn.resize(A.nonZeros());
        jcn.resize(A.nonZeros());
        a.resize(A.nonZeros());
    }
    else
    {
        irn.resize(A.nonZeros());
        jcn.resize(A.nonZeros());
    }

    // Decompress the sparse matrix
    for (auto k = 0, nonzeros = 0; k < A.outerSize(); ++k)
    {
        for (SparseMatrix::InnerIterator it(A, k); it; ++it, ++nonzeros)
        {
            a[nonzeros] = it.value();
            irn[nonzeros] = it.row() + 1;
            jcn[nonzeros] = it.col() + 1;
        }
    }
}

void MUMPSLU::solve(SparseMatrix const& A, Vector& x, Vector const& b)
{
    auto start = std::chrono::high_resolution_clock::now();

    x = b;

    using MUMPSAdapter = MUMPSWrapper<double>;

    MUMPSAdapter::MUMPS_STRUC_C info;

    info.job = Job::Initialization;

    // Par determines parallel data location
    // 0  : Host is not involved with computes
    // 1  : Host is involved
    info.par = 1;
    info.comm_fortran = -987654;

    // Symmetric matrix flag
    info.sym = MatrixProperty::Unsymmetric;

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
    info.nz = A.nonZeros();

    // Convert to MUMPS representation of the matrix

    // Temporary matrix storage
    irn.resize(A.nonZeros());
    jcn.resize(A.nonZeros());
    a.resize(A.nonZeros());

    // Decompress the sparse matrix
    for (auto k = 0, nonzeros = 0; k < A.outerSize(); ++k)
    {
        for (SparseMatrix::InnerIterator it(A, k); it; ++it, ++nonzeros)
        {
            a[nonzeros] = it.value();
            irn[nonzeros] = it.row() + 1;
            jcn[nonzeros] = it.col() + 1;
        }
    }

    info.a = const_cast<double*>(A.coeffs().data());

    info.irn = irn.data();
    info.jcn = jcn.data();

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

void MUMPSLLT::solve(SparseMatrix const& A, Vector& x, Vector const& b)
{
    auto start = std::chrono::high_resolution_clock::now();

    x = b;

    using MUMPSAdapter = MUMPSWrapper<double>;

    MUMPSAdapter::MUMPS_STRUC_C info;

    info.job = Job::Initialization;

    // Par determines parallel data location
    // 0  : Host is not involved with computes
    // 1  : Host is involved
    info.par = 1;
    info.comm_fortran = -987654;

    // Symmetric matrix flag
    info.sym = MatrixProperty::SPD;

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
    info.nz = A.nonZeros();

    // Convert to MUMPS representation of the matrix

    // Temporary matrix storage
    irn.resize(A.nonZeros());
    jcn.resize(A.nonZeros());
    a.resize(A.nonZeros());

    // Decompress the sparse matrix
    for (auto k = 0, nonzeros = 0; k < A.outerSize(); ++k)
    {
        for (SparseMatrix::InnerIterator it(A, k); it; ++it, ++nonzeros)
        {
            a[nonzeros] = it.value();
            irn[nonzeros] = it.row() + 1;
            jcn[nonzeros] = it.col() + 1;
        }
    }

    info.a = const_cast<double*>(A.coeffs().data());

    info.irn = irn.data();
    info.jcn = jcn.data();

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
}
