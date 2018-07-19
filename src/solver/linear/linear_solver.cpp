
// Parallelisation of the Conjugate Gradient algorithm is good
#define NEON_PARALLEL_EIGEN_SOLVERS

#include "linear_solver.hpp"

#include "exceptions.hpp"
#include "simulation_parser.hpp"

#ifdef ENABLE_OPENMP
#include <omp.h>
#endif

#include <cfenv>
#include <chrono>
#include <iostream>

namespace neon
{
iterative_linear_solver::iterative_linear_solver(double const residual_tolerance)
    : residual_tolerance{residual_tolerance}
{
}

iterative_linear_solver::iterative_linear_solver(std::int32_t const max_iterations)
    : max_iterations{max_iterations}
{
}

iterative_linear_solver::iterative_linear_solver(double const residual_tolerance,
                                                 std::int32_t const max_iterations)
    : residual_tolerance{residual_tolerance}, max_iterations(max_iterations)
{
}

void conjugate_gradient::solve(sparse_matrix const& A, vector& x, vector const& b)
{
#ifdef ENABLE_OPENMP
    omp_set_num_threads(simulation_parser::threads);
#endif

    std::feclearexcept(FE_ALL_EXCEPT);

    Eigen::ConjugateGradient<sparse_matrix, Eigen::Lower | Eigen::Upper> pcg;

    pcg.setTolerance(residual_tolerance);
    pcg.setMaxIterations(max_iterations);

    pcg.compute(A);

    x = pcg.solve(b);

    std::cout << std::string(6, ' ') << "Conjugate Gradient iterations: " << pcg.iterations()
              << " (max. " << max_iterations << "), estimated error: " << pcg.error() << " (min. "
              << residual_tolerance << ")\n";

    if (std::fetestexcept(FE_INVALID))
    {
        throw computational_error("Floating point error reported\n");
    }

    if (pcg.iterations() >= max_iterations)
    {
        throw computational_error("Conjugate gradient solver maximum iterations "
                                  "reached");
    }
}

void biconjugate_gradient_stabilised::solve(sparse_matrix const& A, vector& x, vector const& b)
{
    std::feclearexcept(FE_ALL_EXCEPT);

    Eigen::BiCGSTAB<sparse_matrix> bicgstab; //, Eigen::IncompleteLUT<double>

    bicgstab.setTolerance(residual_tolerance);
    bicgstab.setMaxIterations(max_iterations);

    bicgstab.compute(A);

    x = bicgstab.solve(b);

    if (std::fetestexcept(FE_INVALID))
    {
        throw computational_error("Floating point error reported\n");
    }

    if (bicgstab.iterations() >= max_iterations)
    {
        throw computational_error("Conjugate gradient solver maximum iterations "
                                  "reached\n");
    }

    std::cout << std::string(6, ' ') << "Conjugate Gradient iterations: " << bicgstab.iterations()
              << " (max. " << max_iterations << "), estimated error: " << bicgstab.error()
              << " (min. " << residual_tolerance << ")\n";
}

void SparseLU::solve(sparse_matrix const& A, vector& x, vector const& b)
{
    if (build_sparsity_pattern)
    {
        lu.analyzePattern(A);
        build_sparsity_pattern = false;
    }
    lu.factorize(A);
    x = lu.solve(b);
}

void SparseLLT::solve(sparse_matrix const& A, vector& x, vector const& b)
{
    if (build_sparsity_pattern)
    {
        llt.analyzePattern(A);
        build_sparsity_pattern = false;
    }
    llt.factorize(A);
    x = llt.solve(b);
}
}
