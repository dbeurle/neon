
// Parallelisation of the Conjugate Gradient algorithm is good
#define NEON_PARALLEL_EIGEN_SOLVERS

#include "LinearSolver.hpp"

#include <omp.h>

// Must be here
#include <Eigen/PaStiXSupport>

#include <Eigen/Sparse>

// Mumps includes
#include <MUMPS/dmumps_c.h>
#include <MUMPS/smumps_c.h>

#include <chrono>
#include <iostream>
#include <thread>

namespace neon
{
template <typename T>
struct MUMPSWrapper
{
};

template <>
struct MUMPSWrapper<float>
{
    typedef SMUMPS_STRUC_C MUMPS_STRUC_C;
    static void mumps_c(MUMPS_STRUC_C& info) { smumps_c(&info); }
};

template <>
struct MUMPSWrapper<double>
{
    typedef DMUMPS_STRUC_C MUMPS_STRUC_C;
    static void mumps_c(MUMPS_STRUC_C& info) { dmumps_c(&info); }
};

LinearSolver::LinearSolver() : solverParam(1.0e-5, 1000) {}

void PaStiX::solve(const SparseMatrix& A, Vector& x, const Vector& b)
{
    auto start = std::chrono::high_resolution_clock::now();

    Eigen::PastixLLT<Eigen::SparseMatrix<double>, Eigen::Upper> pastix;

    // Verbosity
    pastix.iparm(3) = 1;

    // Number of threads
    pastix.iparm(34) = std::thread::hardware_concurrency();

    // Number of Cuda devices
    // pastix.iparm(64) = 1;

    pastix.compute(A);

    x = pastix.solve(b);

    std::cout << std::string(6, ' ') << "Analysis step " << pastix.dparm(18) << " seconds\n";
    std::cout << std::string(6, ' ') << "Predicted factorisation time " << pastix.dparm(19) << "s\n";
    std::cout << std::string(6, ' ') << "Factorisation " << pastix.dparm(20) << "s\n";
    std::cout << std::string(6, ' ') << "Time for solve " << pastix.dparm(21) << "s\n";
    std::cout << std::string(6, ' ') << "GigaFLOPS during factorisation "
              << pastix.dparm(22) / 1.0e9 << "\n";
    std::cout << std::string(6, ' ') << "MegaFLOPS during solve " << pastix.dparm(23) / 1.0e6 << "\n";

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;
}

void MUMPS::solve(SparseMatrix const& A, Vector& x, Vector const& b)
{
    auto start = std::chrono::high_resolution_clock::now();

    x = b;

    // Check for a square matrix
    assert(A.rows() == A.cols() == x.size() == b.size());

    MUMPSWrapper<double>::MUMPS_STRUC_C info;

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

    MUMPSWrapper<double>::mumps_c(info);

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
    MUMPSWrapper<double>::mumps_c(info);

    if (info.info[0] < 0) throw std::runtime_error("Error in analysis phase of MUMPS solver\n");

    // Factorization phase
    info.job = 2;
    MUMPSWrapper<double>::mumps_c(info);

    if (info.info[0] < 0) throw std::runtime_error("Error in factorisation phase\n");

    info.rhs = x.data();
    info.nrhs = 1;
    info.lrhs = info.n;

    // Back substitution
    info.job = 3;
    MUMPSWrapper<double>::mumps_c(info);

    if (info.info[0] < 0)
        throw std::runtime_error("Error in back-substitution phase of MUMPS solver\n");

    // Take out the trash
    info.job = -2;
    MUMPSWrapper<double>::mumps_c(info);

    if (info.info[0] < 0) throw std::runtime_error("Error in cleanup phase\n");

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;

    std::cout << "    Linear solver took " << elapsed_seconds.count() << "s\n";
}

void SparseLU::solve(SparseMatrix const& A, Vector& x, Vector const& b)
{
    Eigen::SparseLU<SparseMatrix, Eigen::AMDOrdering<int>> sparseLU;
    sparseLU.compute(A);
    x = sparseLU.solve(b);
}

pCG::pCG(double tol) { solverParam.tolerance = tol; }

pCG::pCG(int maxIter) { solverParam.max_iterations = maxIter; }

pCG::pCG(double tol, int maxIter)
{
    solverParam.max_iterations = maxIter;
    solverParam.tolerance = tol;
}

void pCG::solve(SparseMatrix const& A, Vector& x, const Vector& b)
{
    omp_set_num_threads(std::thread::hardware_concurrency());
    // omp_set_num_threads(1);

    Eigen::ConjugateGradient<SparseMatrix, Eigen::Lower | Eigen::Upper> pcg;

    pcg.setTolerance(LinearSolver::solverParam.tolerance);
    pcg.setMaxIterations(LinearSolver::solverParam.max_iterations);

    pcg.compute(A);

    x = pcg.solveWithGuess(b, x);

    std::cout << std::string(6, ' ') << "Conjugate Gradient iterations: " << pcg.iterations()
              << " (max. " << solverParam.max_iterations << "), estimated error: " << pcg.error()
              << " (min. " << solverParam.tolerance << ")\n";

    if (pcg.iterations() >= solverParam.max_iterations)
    {
        throw std::runtime_error("Conjugate gradient solver maximum iterations reached.  Try "
                                 "increasing the maximum number of iterations or use a different "
                                 "solver\n");
    }
    // std::cout << "    Linear solver took " << elapsed_seconds.count() << "s\n";
}

BiCGSTAB::BiCGSTAB(double tol) { solverParam.tolerance = tol; }

BiCGSTAB::BiCGSTAB(int maxIter) { solverParam.max_iterations = maxIter; }

BiCGSTAB::BiCGSTAB(double tol, int maxIter)
{
    solverParam.max_iterations = maxIter;
    solverParam.tolerance = tol;
}

void BiCGSTAB::solve(const SparseMatrix& A, Vector& x, const Vector& b)
{
    Eigen::BiCGSTAB<SparseMatrix> bicgstab; //, Eigen::IncompleteLUT<double>

    bicgstab.setTolerance(solverParam.tolerance);
    bicgstab.setMaxIterations(solverParam.max_iterations);

    auto start = std::chrono::high_resolution_clock::now();

    bicgstab.compute(A);

    x = bicgstab.solve(b);

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;

    std::cout << "#iterations:     " << bicgstab.iterations();
    std::cout << "estimated error: " << bicgstab.error();
}
}
