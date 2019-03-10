
// Parallelisation of the Conjugate Gradient algorithm is good
#define NEON_PARALLEL_EIGEN_SOLVERS

#include "linear_solver.hpp"

#include "exceptions.hpp"
#include "simulation_parser.hpp"
#include "graph/cuthill_mckee.hpp"
#include "graph/bandwidth.hpp"

#ifdef ENABLE_OPENMP
#include <omp.h>
#endif

#include <cfenv>
#include <chrono>
#include <iostream>
#include <fstream>

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

    if (build_sparsity_pattern)
    {
        auto const start = std::chrono::steady_clock::now();

        reverse_cuthill_mcgee reordering(A);

        reordering.compute();

        auto const& permutation = reordering.permutation();

        P.indices().resize(permutation.size());

        std::copy(begin(permutation), end(permutation), P.indices().data());

        build_sparsity_pattern = false;

        std::chrono::duration<double> const elapsed_seconds = std::chrono::steady_clock::now()
                                                              - start;

        std::cout << std::string(6, ' ') << "Reordering took " << elapsed_seconds.count() << "s\n";
    }

    auto const p_start = std::chrono::steady_clock::now();

    permuted_matrix = P.transpose() * A * P;
    permuted_rhs = P.transpose() * b;

    std::cout << std::string(8, ' ') << "Bandwidth reduced from " << compute_bandwidth(A) << " to "
              << compute_bandwidth(permuted_matrix) << "\n";

    std::chrono::duration<double> const p_time = std::chrono::steady_clock::now() - p_start;

    std::cout << std::string(8, ' ') << "Permutation took " << p_time.count() << "s\n";

    Eigen::ConjugateGradient<sparse_matrix, Eigen::Lower | Eigen::Upper> pcg;

    pcg.setTolerance(residual_tolerance);
    pcg.setMaxIterations(max_iterations);

    auto const start = std::chrono::steady_clock::now();

    x = P * pcg.compute(permuted_matrix).solve(permuted_rhs);

    auto const end = std::chrono::steady_clock::now();
    std::chrono::duration<double> const elapsed_seconds = end - start;

    std::cout << std::string(6, ' ') << "Conjugate gradient took " << elapsed_seconds.count()
              << "s, iterations: " << pcg.iterations() << " (max. " << max_iterations
              << "), estimated error: " << pcg.error() << " (min. " << residual_tolerance << ")\n";

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

    std::cout << std::string(6, ' ') << "Conjugate gradient iterations: " << bicgstab.iterations()
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
