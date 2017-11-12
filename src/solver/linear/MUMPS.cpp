
#include "MUMPS.hpp"
#include "Exceptions.hpp"

#include <Eigen/Sparse>

#include <chrono>
#include <iostream>

namespace neon
{
MUMPS::MUMPS(MatrixProperty const symmetric_flag)
{
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
}

MUMPS::~MUMPS()
{
    // Take out the trash
    info.job = Job::Terminate;
    MUMPSAdapter::mumps_c(info);
}

void MUMPS::internal_solve(SparseMatrix const& A, Vector& x, Vector const& b)
{
    auto start = std::chrono::high_resolution_clock::now();

    x = b;

    info.n = A.rows();
    info.nz = coefficients.size();

    info.a = coefficients.data();
    info.irn = rows.data();
    info.jcn = cols.data();

    if (build_sparsity_pattern)
    {
        // Analysis phase
        info.job = Job::Analysis;
        MUMPSAdapter::mumps_c(info);

        if (info.info[0] < 0)
        {
            throw computational_error("Error in analysis phase of MUMPS solver\n");
        }
        build_sparsity_pattern = false;
    }

    // Factorization phase
    info.job = Job::Factorisation;
    MUMPSAdapter::mumps_c(info);

    if (info.info[0] < 0)
    {
        throw computational_error("Error in factorisation phase of MUMPS solver\n");
    }

    info.rhs = x.data();
    info.nrhs = 1;
    info.lrhs = info.n;

    info.job = Job::BackSubstitution;
    MUMPSAdapter::mumps_c(info);

    if (info.info[0] < 0)
    {
        throw computational_error("Error in back substitution phase of MUMPS solver\n");
    }

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;

    std::cout << std::string(6, ' ') << "MUMPS solver took " << elapsed_seconds.count() << "s\n";
}

void MUMPSLLT::allocate_coordinate_format_storage(SparseMatrix const& A)
{
    coefficients.clear();
    coefficients.reserve(A.nonZeros());

    if (build_sparsity_pattern)
    {
        rows.reserve(A.nonZeros());
        cols.reserve(A.nonZeros());

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
        rows.shrink_to_fit();
        cols.shrink_to_fit();
    }
    else
    {
        // Only update the non-zero numerical values
        for (auto k = 0; k < A.outerSize(); ++k)
        {
            for (SparseMatrix::InnerIterator it(A, k); it; ++it)
            {
                if (it.col() >= it.row())
                {
                    coefficients.emplace_back(it.value());
                }
            }
        }
    }
}

void MUMPSLLT::solve(SparseMatrix const& A, Vector& x, Vector const& b)
{
    this->allocate_coordinate_format_storage(A);
    internal_solve(A, x, b);
}

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

void MUMPSLU::solve(SparseMatrix const& A, Vector& x, Vector const& b)
{
    this->allocate_coordinate_format_storage(A);
    internal_solve(A, x, b);
}
}
