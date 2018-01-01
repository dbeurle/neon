
#include "distributed_pastix.hpp"

#include "SimulationControl.hpp"

#include <chrono>
#include <termcolor/termcolor.hpp>

namespace neon
{
distributed_pastix_ldlt::distributed_pastix_ldlt()
{
    // Verbosity
    // ldlt.iparm(3) = 0;

    // Number of threads
    iparm[IPARM_THREAD_NBR] = SimulationControl::threads;

    // Number of Cuda devices
    // ldlt.iparm(64) = 1;
}

distributed_pastix_ldlt::solve(SparseMatrix const& A, vector& x, vector const& b)
{
    x = b;

    if (build_sparsity_pattern)
    {
        // Permutation tabular
        dpastix(&pastix_data,
                MPI_COMM_WORLD,
                A.rows(),
                colptr.data(),
                rows.data(),
                A.valuePtr(),
                mapping.data(),
                perm.data(),
                invp.data(),
                x.data(),
                1,
                iparm.data(),
                dparm.data());

        build_sparsity_pattern = false;
    }

    // Customise parameters
    iparm[IPARM_SYM] = API_SYM_YES;
    iparm[IPARM_FACTORIZATION] = API_FACT_LDLT;
    iparm[IPARM_MATRIX_VERIFICATION] = API_NO;
    iparm[IPARM_VERBOSE] = 2;
    iparm[IPARM_ORDERING] = API_ORDER_SCOTCH;
    iparm[IPARM_RHS_MAKING] = API_RHS_B;

    // Re-read parameters to set IPARM/DPARM
    if (EXIT_FAILURE == get_idparm(argc, argv, iparm.data(), dparm.data()))
    {
        return EXIT_FAILURE;
    }

    iparm[IPARM_START_TASK] = API_TASK_ORDERING;
    iparm[IPARM_END_TASK] = API_TASK_CLEAN;

    // Copy the right hand side as it will be replaced by solution
    auto rhssaved = rhs;

    pastix_int_t number_of_dofs = mpi::all_reduce(local_number_of_dofs, mpi::sum{});

    // right hand side (save, global)
    std::vector<pastix_float_t> rhssaved_g(number_of_dofs, 0.0);
    for (auto i = 0; i < local_number_of_dofs; i++)
    {
        rhssaved_g[mapping[i] - 1] = rhssaved[i];
    }

    rhssaved_g = mpi::all_reduce(rhssaved_g, mpi::sum{});

    // Call pastix
    dpastix(&pastix_data,
            MPI_COMM_WORLD,
            local_number_of_dofs,
            colptr.data(),
            rows.data(),
            values.data(),
            mapping.data(),
            perm.data(),
            nullptr, // No need to allocate invp in dpastix
            rhs.data(),
            1,
            iparm.data(),
            dparm.data());

    for (auto x : rhs)
    {
        if (std::abs(x - 1.0) > 0.001)
        {
            std::cout << x << std::endl;
            std::cout << "Solution not correct!" << std::endl;
            mpi::finalise();
        }
    }
}

void distributed_pastix_ldlt::analyse_pattern() {}
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
