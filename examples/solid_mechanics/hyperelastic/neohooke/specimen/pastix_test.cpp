
#include <array>
#include <vector>
#include <iostream>

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <unsupported/Eigen/SparseExtra>

extern "C" {
#include <pastix.h>
}

int main()
{
    using llt_solver = Eigen::SimplicialLLT<Eigen::SparseMatrix<double, Eigen::ColMajor>, Eigen::Lower>;

    using sparse_matrix = Eigen::SparseMatrix<double, Eigen::RowMajor>;
    using vector = Eigen::VectorXd;

    sparse_matrix A;
    Eigen::loadMarket(A, "sparse_matrix.mtx");
    A.finalize();

    vector b;
    Eigen::loadMarketVector(b, "right_hand_side.mtx");

    // Solve using Eigen built in solvers
    // vector x = llt_solver{}.compute(A).solve(b);
    // std::cout << "Eigen solution:\n" << (A * x - b).norm() << std::endl;

    vector x = b;

    assert(x.rows() == A.cols());
    assert(x.rows() == b.rows());
    assert(x.cols() == 1);
    assert(b.cols() == 1);

    constexpr auto nrhs = 1;

    // Integer in/out parameters for pastix
    std::array<pastix_int_t, IPARM_SIZE> iparm;
    // Floating in/out parameters for pastix
    std::array<double, DPARM_SIZE> dparm;

    // Initialize parameters to default values
    pastixInitParam(iparm.data(), dparm.data());

    std::cout << "Integer parameters\n";
    for (auto const& i : iparm)
    {
        std::cout << i << "\n";
    }

    std::cout << "Double parameters\n";
    for (auto const& d : dparm)
    {
        std::cout << d << "\n";
    }

    auto spm = spmatrix_t{};

    spm.mtxtype = SpmGeneral;
    spm.flttype = SpmDouble;
    spm.fmttype = A.IsRowMajor ? SpmCSR : SpmCSC;

    spm.gN = 0;
    spm.n = A.rows();
    spm.gnnz = 0;
    spm.nnz = A.nonZeros();

    spm.gNexp = 0;
    spm.nexp = 0;
    spm.gnnzexp = 0;
    spm.nnzexp = 0;

    spm.dof = 1;
    spm.dofs = nullptr;
    spm.layout = A.IsRowMajor ? SpmRowMajor : SpmColMajor;

    spm.rowptr = const_cast<pastix_int_t*>(A.IsRowMajor ? A.innerIndexPtr() : A.outerIndexPtr());
    spm.colptr = const_cast<pastix_int_t*>(A.IsRowMajor ? A.outerIndexPtr() : A.innerIndexPtr());
    spm.loc2glob = nullptr;
    spm.values = const_cast<double*>(A.valuePtr());

    spmUpdateComputedFields(&spm);

    spmPrintInfo(&spm, stdout);

    iparm[IPARM_FACTORIZATION] = PastixFactLDLT;

    // Pointer to the storage structure required by pastix
    pastix_data_t* pastix_data = nullptr;

    // Startup PaStiX
    pastixInit(&pastix_data, MPI_COMM_WORLD, iparm.data(), dparm.data());

    // Perform ordering, symbolic factorization, and analyze steps
    pastix_task_analyze(pastix_data, &spm);

    // Perform the numerical factorisation
    pastix_task_numfact(pastix_data, &spm);

    x = b;

    // Solve the linear system
    pastix_task_solve(pastix_data, nrhs, const_cast<double*>(x.data()), spm.n);

    // Iteratively refine the solution
    pastix_task_refine(pastix_data,
                       spm.n,
                       nrhs,
                       const_cast<double*>(b.data()),
                       b.rows(),
                       x.data(),
                       x.rows());

    pastixFinalize(&pastix_data);

    std::cout << "Solution norm : " << (A * x - b).norm() << "\n";

    return 0;
}
