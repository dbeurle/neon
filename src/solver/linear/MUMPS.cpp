
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

void MUMPS::solve(SparseMatrix const& A, Vector& x, Vector const& b)
{
    auto start = std::chrono::high_resolution_clock::now();

    x = b;

    // Check for a square matrix
    assert(A.rows() == A.cols());
    assert(x.size() == b.size());
    assert(A.cols() == x.size());

    using MUMPSadapter = MUMPSWrapper<double>;

    MUMPSadapter::MUMPS_STRUC_C info;

    // Jobs in MUMPS use the following:
    // -1  Solver initialization
    // -2  Termination
    // 1   Perform the analysis
    // 2   Perform the factorization
    // 3   Compute the solution
    // 4   Job = 1 && Job = 2
    // 5   Job = 2 && Job = 3
    // 6   Job = 1 && Job = 2 && Job = 3
    info.job = -1;

    // Par determines parallel data location
    // 0  : Host is not involved with computes
    // 1  : Host is involved
    info.par = 1;
    info.comm_fortran = -987654;

    // Symmetric matrix flag
    // 0  A is unsymmetric
    // 1  A is SPD
    // 2  A is general symmetric
    info.sym = 0;

    MUMPSadapter::mumps_c(info);

    // Verbosity levels
    info.icntl[0] = 1;
    info.icntl[1] = 1;
    info.icntl[2] = 1;
    info.icntl[3] = 1;

    // Ordering algorithm
    // 0 	AMD
    // 2 	AMF
    // 3 	Scotch
    // 4	Pord
    // 5	Metis
    // 6    QAMD
    // 7    Automatic based on matrix
    info.icntl[6] = 7;

    // Iterative refinement
    info.icntl[9] = 100;

    // Compute the residual
    // 0	No residuals
    // 1 	Expensive (condition number etc)
    // 2 	Cheap residuals
    info.icntl[10] = 2;

    info.n = A.rows();
    info.nz = A.nonZeros();

    // Convert to MUMPS representation of the matrix

    // Temporary matrix storage
    std::vector<int> irn(info.nz);
    std::vector<int> jcn(info.nz);
    std::vector<double> a(info.nz);

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

    info.a = a.data();
    info.irn = irn.data();
    info.jcn = jcn.data();

    // Analysis phase
    info.job = 1;
    MUMPSadapter::mumps_c(info);

    if (info.info[0] < 0)
        throw computational_error("Error in analysis phase of MUMPS solver\n");

    // Factorization phase
    info.job = 2;
    MUMPSadapter::mumps_c(info);

    if (info.info[0] < 0) throw computational_error("Error in factorisation phase\n");

    info.rhs = x.data();
    info.nrhs = 1;
    info.lrhs = info.n;

    // Back substitution
    info.job = 3;
    MUMPSadapter::mumps_c(info);

    if (info.info[0] < 0)
        throw computational_error("Error in back-substitution phase of MUMPS solver\n");

    // Take out the trash
    info.job = -2;
    MUMPSadapter::mumps_c(info);

    if (info.info[0] < 0) throw computational_error("Error in cleanup phase\n");

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;

    std::cout << std::string(6, ' ') << "MUMPS solver took " << elapsed_seconds.count()
              << "s\n";
}
}
