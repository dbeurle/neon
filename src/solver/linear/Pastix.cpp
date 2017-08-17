
#include "Pastix.hpp"

#include "SimulationControl.hpp"

#include <chrono>
#include <termcolor/termcolor.hpp>

#include <Eigen/PaStiXSupport>

namespace neon
{
void PaStiX::solve(SparseMatrix const& A, Vector& x, Vector const& b)
{
    auto start = std::chrono::high_resolution_clock::now();

    Eigen::PastixLLT<Eigen::SparseMatrix<double>, Eigen::Upper> pastix;

    // Verbosity
    pastix.iparm(3) = 0;

    // Number of threads
    pastix.iparm(34) = SimulationControl::threads;

    // Number of Cuda devices
    // pastix.iparm(64) = 1;

    pastix.compute(A);

    x = pastix.solve(b);

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;
    std::cout << std::string(6, ' ') << "PaStiX LLT direct solver took " << elapsed_seconds.count()
              << "s\n";
}
}

// FIXME When new pastix comes out, test this against it.  Currently segfaults
// upon use of A after passing the col and row pointer to PaStiX

//     extern "C" {
// #include <pastix.h>
// #include <pastix/nompi.h>
// #include <spm.h>
//     }

// SparseMatrix A = A_in;
//
// constexpr int nrhs = 1;
//
// std::array<pastix_int_t, IPARM_SIZE> iparm; // Integer in/out parameters for pastix
// std::array<double, DPARM_SIZE> dparm;       // Floating in/out parameters for
// pastix
//
// // Initialize parameters to default values
// pastixInitParam(iparm.data(), dparm.data());
//
// std::cout << "Integer parameters\n";
// for (auto const& i : iparm) std::cout << i << std::endl;
//
// std::cout << "Double parameters\n";
// for (auto const& d : dparm) std::cout << d << std::endl;
//
// auto* spm = spmNew(PastixSymmetric,
//                    PastixDouble,
//                    A.IsRowMajor ? PastixCSR : PastixCSC,
//                    A.rows(),
//                    A.nonZeros(),
//                    const_cast<pastix_int_t*>(A.IsRowMajor ? A.outerIndexPtr()
//                                                           : A.innerIndexPtr()),
//                    const_cast<pastix_int_t*>(A.IsRowMajor ? A.innerIndexPtr()
//                                                           : A.outerIndexPtr()),
//                    const_cast<double*>(A.valuePtr()),
//                    nullptr,
//                    1,
//                    A.IsRowMajor ? PastixRowMajor : PastixColMajor,
//                    nullptr);
//
// spmPrintInfo(spm, stdout);
//
// iparm[IPARM_FACTORIZATION] = PastixFactLLT;
//
// // Pointer to the storage structure required by pastix
// pastix_data_t* pastix_data = nullptr;
//
// // Startup PaStiX
// pastixInit(&pastix_data, MPI_COMM_WORLD, iparm.data(), dparm.data());
//
// // Perform ordering, symbolic factorization, and analyze steps
// pastix_task_analyze(pastix_data, spm);
//
// // Perform the numerical factorisation
// pastix_task_numfact(pastix_data, spm);
//
// x = b;
//
// // Solve the linear system
// pastix_task_solve(pastix_data, nrhs, const_cast<double*>(x.data()), spm->n);
//
// // Iteratively refine the solution
// // pastix_task_refine(pastix_data, x.data(), nrhs, const_cast<double*>(b.data()));
//
// spmFree(spm);
//
// pastixFinalize(&pastix_data);
//
// std::cout << "Pastix finalized...!" << std::endl;
//
// std::cout << "Solution norm : " << (A_in * x - b).norm() << std::endl;
//
// for (int k = 0; k < A_in.outerSize(); ++k)
//     for (SparseMatrix::InnerIterator it(A_in, k); it; ++it)
//     {
//         // it.value();
//         std::cout << "(" << it.row() << ", " << it.col() << ") with index "
//                   << it.index() << "\n"; // row index
//     }
//
// std::cout << A_in << std::endl;
