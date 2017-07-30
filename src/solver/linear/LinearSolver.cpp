
// Parallelisation of the Conjugate Gradient algorithm is good
#define NEON_PARALLEL_EIGEN_SOLVERS

#include "LinearSolver.hpp"

#include "SimulationControl.hpp"

#ifdef ENABLE_OPENMP
#include <omp.h>
#endif

#include <Eigen/Sparse>

#include <chrono>
#include <iostream>
#include <thread>

namespace neon
{
LinearSolver::LinearSolver() : solverParam(1.0e-5, 1000) {}

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
#ifdef ENABLE_OPENMP
    omp_set_num_threads(SimulationControl::threads);
#endif

    Eigen::ConjugateGradient<SparseMatrix, Eigen::Lower | Eigen::Upper> pcg;

    pcg.setTolerance(LinearSolver::solverParam.tolerance);
    pcg.setMaxIterations(LinearSolver::solverParam.max_iterations);

    pcg.compute(A);

    x = pcg.solveWithGuess(b, x);

    std::cout << std::string(6, ' ')
              << "Conjugate Gradient iterations: " << pcg.iterations() << " (max. "
              << solverParam.max_iterations << "), estimated error: " << pcg.error()
              << " (min. " << solverParam.tolerance << ")\n";

    if (pcg.iterations() >= solverParam.max_iterations)
    {
        throw std::runtime_error(
            "Conjugate gradient solver maximum iterations reached.  Try "
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
