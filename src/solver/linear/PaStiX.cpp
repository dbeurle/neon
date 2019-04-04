
#include "PaStiX.hpp"

#include "simulation_parser.hpp"

#include <chrono>
#include <termcolor/termcolor.hpp>

namespace neon
{
PaStiXLDLT::PaStiXLDLT()
{
    // Verbosity
    ldlt.iparm(3) = 0;

    // Number of threads
    ldlt.iparm(34) = simulation_parser::threads;

    // Number of Cuda devices
    // ldlt.iparm(64) = 1;
}

void PaStiXLDLT::solve(sparse_matrix const& A, vector& x, vector const& b)
{
    auto const start = std::chrono::steady_clock::now();

    if (build_sparsity_pattern)
    {
        ldlt.analyzePattern(A);
        build_sparsity_pattern = false;
    }

    ldlt.factorize(A);

    x = ldlt.solve(b);

    auto const end = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;
    std::cout << std::string(6, ' ') << "PaStiX LDLT direct solver took " << elapsed_seconds.count()
              << "s\n";
}

PaStiXLU::PaStiXLU()
{
    // Verbosity
    lu.iparm(3) = 0;

    // Number of threads
    lu.iparm(34) = simulation_parser::threads;

    // Number of Cuda devices
    // ldlt.iparm(64) = 1;
}

void PaStiXLU::solve(sparse_matrix const& A, vector& x, vector const& b)
{
    auto start = std::chrono::steady_clock::now();

    if (build_sparsity_pattern)
    {
        lu.analyzePattern(A);
        build_sparsity_pattern = false;
    }

    lu.factorize(A);

    x = lu.solve(b);

    auto end = std::chrono::steady_clock::now();
    std::chrono::duration<double> const elapsed_seconds = end - start;
    std::cout << std::string(6, ' ') << "PaStiX LU direct solver took " << elapsed_seconds.count()
              << "s\n";
}
}
