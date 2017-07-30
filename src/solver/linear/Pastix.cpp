
#include "Pastix.hpp"

#include "SimulationControl.hpp"

// Must be here
#include <Eigen/PaStiXSupport>

#include <chrono>
#include <thread>

namespace neon
{
void PaStiX::solve(const SparseMatrix& A, Vector& x, const Vector& b)
{
    auto start = std::chrono::high_resolution_clock::now();

    Eigen::PastixLLT<Eigen::SparseMatrix<double>, Eigen::Upper> pastix;

    // Verbosity
    pastix.iparm(3) = 1;

    // Number of threads
    pastix.iparm(34) = SimulationControl::threads;

    // Number of Cuda devices
    // pastix.iparm(64) = 1;

    pastix.compute(A);

    x = pastix.solve(b);

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;
}
}
