
#include "LinearSolver.hpp"

// Must be here
#include <Eigen/PaStiXSupport>

#include <Eigen/Sparse>

// Mumps includes
#include <MUMPS/dmumps_c.h>
#include <MUMPS/smumps_c.h>

#include <boost/timer/timer.hpp>
#include <iostream>

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

LinearSolver::LinearSolver() : solverParam(1e-6, 1000)
{
    solverParam.max_iterations = 1000;
    solverParam.tolerance = 1e-5;
}

void PaStiX::solve(const SparseMatrix& A, Vector& x, const Vector& b)
{
    boost::timer::cpu_timer timer;

    Eigen::PastixLU<Eigen::SparseMatrix<double>> pastix;

    // Number of threads
    pastix.iparm(34) = 2;

    // Number of Cuda devices
    // pastix.iparm(64) = 1;

    pastix.compute(A);

    x = pastix.solve(b);

    std::cout << "\tAnalysis step " << pastix.dparm(18) << " seconds\n";
    std::cout << "\tPredicted factorization time " << pastix.dparm(19) << " seconds\n";
    std::cout << "\tFactorization " << pastix.dparm(20) << " seconds\n";
    std::cout << "\tTime for solve " << pastix.dparm(21) << " seconds\n";
    std::cout << "\tGigaFLOPS during factorization " << pastix.dparm(22) / 1e9 << "\n";
    std::cout << "\tMegaFLOPS during solve " << pastix.dparm(23) / 1e6 << "\n";
    std::cout << "Solution time: " << timer.format() << "\n";
}

void MUMPS::solve(const SparseMatrix& A, Vector& x, const Vector& b)
{
    boost::timer::cpu_timer timer;

    x = b;

    auto const rowcols = std::make_pair(A.rows(), A.cols());

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

    info.n = rowcols.first;
    info.nz = A.nonZeros();

    // Convert to MUMPS representation of the matrix

    // Temporary matrix storage
    int* irn = new int[info.nz];
    int* jcn = new int[info.nz];
    double* a = new double[info.nz];

    // Decompress the sparse matrix
    for (auto k = 0, nonzeros = 0; k < A.outerSize(); ++k)
    {
        for (SparseMatrix::InnerIterator it(A, k); it; ++it)
        {
            a[nonzeros] = it.value();
            irn[nonzeros] = it.row() + 1;
            jcn[nonzeros] = it.col() + 1;
            ++nonzeros;
        }
    }

    info.a = a;
    info.irn = irn;
    info.jcn = jcn;

    // Analysis phase
    info.job = 1;
    MUMPSWrapper<double>::mumps_c(info);

    if (info.info[0] < 0)
    {
        std::cout << "Error in analysis phase\n";
        // FIXME throw exception
        std::abort();
    }

    // Factorization phase
    info.job = 2;
    MUMPSWrapper<double>::mumps_c(info);

    if (info.info[0] < 0)
    {
        std::cout << "Error in factorization phase\n";
        // FIXME throw exception
        std::abort();
    }

    info.rhs = x.data();
    info.nrhs = 1;
    info.lrhs = info.n;

    // Back substitution
    info.job = 3;
    MUMPSWrapper<double>::mumps_c(info);
    if (info.info[0] < 0)
    {
        std::cout << "Error in backsubstitution phase\n";
        // FIXME throw exception
        std::abort();
    }

    // Take out the trash
    info.job = -2;
    MUMPSWrapper<double>::mumps_c(info);

    if (info.info[0] < 0)
    {
        std::cout << "Error in cleanup phase\n";
        std::abort();
    }

    std::cout << "Solution time: " << timer.format();

    // Clean up matrix copy
    delete[] a;
    delete[] irn;
    delete[] jcn;
    a = nullptr;
    irn = nullptr;
    jcn = nullptr;
}

void SparseLU::solve(const SparseMatrix& A, Vector& x, const Vector& b)
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

void pCG::solve(const SparseMatrix& A, Vector& x, const Vector& b)
{
    Eigen::ConjugateGradient<SparseMatrix, Eigen::Lower | Eigen::Upper> pcg;

    pcg.setTolerance(LinearSolver::solverParam.tolerance);
    pcg.setMaxIterations(LinearSolver::solverParam.max_iterations);

    std::cout << "Solver tolerance      = " << solverParam.tolerance << "\n";
    std::cout << "Solver max iterations = " << solverParam.max_iterations << "\n";

    boost::timer::cpu_timer timer;

    pcg.compute(A);
    x = pcg.solveWithGuess(b, x);

    std::cout << "Iterations: " << pcg.iterations() << ", estimated error: " << pcg.error() << "\n";
    std::cout << "Solution time: " << timer.format() << "\n";
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

    boost::timer::cpu_timer timer;

    bicgstab.compute(A);

    x = bicgstab.solve(b);

    std::cout << "#iterations:     " << bicgstab.iterations();
    std::cout << "estimated error: " << bicgstab.error();
    std::cout << "Solution time: " << timer.format();
}
}
