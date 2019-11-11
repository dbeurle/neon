
#include "PaStiX.hpp"

#include <chrono>
#include <termcolor/termcolor.hpp>

#include <unsupported/Eigen/SparseExtra>

namespace neon
{
PaStiX::PaStiX() : m_matrix(std::make_unique<spmatrix_t>())
{
    // pastixInitParam(m_integer_parameters.data(), m_float_parameters.data());
    //
    // pastixInit(&m_data, MPI_COMM_WORLD, m_integer_parameters.data(), m_float_parameters.data());
    //
    // m_matrix->mtxtype = SpmGeneral;
    // m_matrix->flttype = SpmDouble;
    // m_matrix->fmttype = SpmCSR;
    //
    // m_matrix->gN = 0;
    // m_matrix->n = 0;
    // m_matrix->gnnz = 0;
    // m_matrix->nnz = 0;
    //
    // m_matrix->gNexp = 0;
    // m_matrix->nexp = 0;
    // m_matrix->gnnzexp = 0;
    // m_matrix->nnzexp = 0;
    //
    // m_matrix->dof = 1;
    // m_matrix->dofs = nullptr;
    // m_matrix->layout = SpmRowMajor;
    //
    // m_matrix->colptr = nullptr;
    // m_matrix->rowptr = nullptr;
    // m_matrix->loc2glob = nullptr;
    // m_matrix->values = nullptr;
}

PaStiX::~PaStiX()
{
    // pastixFinalize(&m_data);
}

void PaStiX::analyse_pattern() { pastix_task_analyze(m_data, m_matrix.get()); }

void PaStiXLDLT::solve(sparse_matrix const& A_in, vector& x_out, vector const& b_in)
{
    auto const start = std::chrono::steady_clock::now();

    sparse_matrix A = A_in;
    vector b = b_in;
    vector x = b;

    Eigen::saveMarket(A, "sparse_matrix.mtx");
    Eigen::saveMarketVector(b, "right_hand_side.mtx");

    // sparse_matrix A;
    // Eigen::loadMarket(A, "sparse_matrix.mtx");
    // A.finalize();
    // vector b;
    // Eigen::loadMarketVector(b, "right_hand_side.mtx");
    // vector x = b;

    assert(A.rows() == A.cols());
    assert(A.isCompressed());

    assert(x.rows() == A.cols());
    assert(x.rows() == b.rows());
    assert(x.cols() == 1);
    assert(b.cols() == 1);

    // Integer in/out parameters for pastix
    std::array<pastix_int_t, IPARM_SIZE> iparm;
    // Floating in/out parameters for pastix
    std::array<double, DPARM_SIZE> dparm;

    // Initialize parameters to default values
    pastixInitParam(iparm.data(), dparm.data());

    // std::cout << "Integer parameters\n";
    // for (auto const& i : iparm)
    // {
    //     std::cout << i << "\n";
    // }
    //
    // std::cout << "Double parameters\n";
    // for (auto const& d : dparm)
    // {
    //     std::cout << d << "\n";
    // }

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

    constexpr auto nrhs = 1;

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
    pastix_task_solve(pastix_data, nrhs, x.data(), spm.n);

    // Iteratively refine the solution
    // pastix_task_refine(pastix_data, spm.n, nrhs, b.data(), b.rows(), x.data(), x.rows());

    pastixFinalize(&pastix_data);

    x_out = x;

    if (std::isnan(b.norm()))
    {
        throw std::runtime_error("Right hand side norm is nan");
    }

    if (std::isnan(x.norm()))
    {
        throw std::runtime_error("Solution norm is nan");
    }

    std::abort();

    // m_matrix->n = A.rows();
    // m_matrix->nnz = A.nonZeros();
    //
    // m_matrix->fmttype = A.IsRowMajor ? SpmCSR : SpmCSC;
    // m_matrix->layout = A.IsRowMajor ? SpmRowMajor : SpmColMajor;
    //
    // m_matrix->colptr = const_cast<pastix_int_t*>(A.IsRowMajor ? A.outerIndexPtr()
    //                                                           : A.innerIndexPtr());
    // m_matrix->rowptr = const_cast<pastix_int_t*>(A.IsRowMajor ? A.innerIndexPtr()
    //                                                           : A.outerIndexPtr());
    // m_matrix->values = const_cast<double*>(A.valuePtr());
    //
    // m_integer_parameters[IPARM_FACTORIZATION] = PastixFactLDLT;
    //
    // if (build_sparsity_pattern)
    // {
    //     spmUpdateComputedFields(m_matrix.get());
    //
    //     spmPrintInfo(m_matrix.get(), stdout);
    //
    //     this->analyse_pattern();
    //     build_sparsity_pattern = false;
    // }
    //
    // auto constexpr nrhs = 1;
    //
    // if (pastix_task_numfact(m_data, m_matrix.get()) != PASTIX_SUCCESS)
    // {
    //     throw std::runtime_error("error in factorisation");
    // }
    //
    // if (pastix_task_solve(m_data, nrhs, x.data(), m_matrix->n) != PASTIX_SUCCESS)
    // {
    //     throw std::runtime_error("error in solve");
    // }
    //
    // if (std::isnan(b.norm()))
    // {
    //     throw std::runtime_error("Right hand side norm is nan");
    // }
    //
    // if (std::isnan(x.norm()))
    // {
    //     throw std::runtime_error("Solution norm is nan");
    // }
    //
    // if (pastix_task_refine(m_data,
    //                        m_matrix->n,
    //                        nrhs,
    //                        const_cast<double*>(b.data()),
    //                        b.rows(),
    //                        x.data(),
    //                        x.rows())
    //     != PASTIX_SUCCESS)
    // {
    //     throw std::runtime_error("error in refine step");
    // }
    //
    // x_out = x;
    //
    // std::cout << "pastix norm in program: " << (A * x - b).norm() << "\n";
    // std::cout << "solution norm: " << x.norm() << "\n";
    //
    // auto const end = std::chrono::steady_clock::now();
    // std::chrono::duration<double> const elapsed_seconds = end - start;
    // std::cout << std::string(6, ' ') << "PaStiX LDLT direct solver took " << elapsed_seconds.count()
    //           << "s\n";
}

PaStiXLU::PaStiXLU() {}

void PaStiXLU::solve(sparse_matrix const& A, vector& x, vector const& b)
{
    auto const start = std::chrono::steady_clock::now();

    assert(A.rows() == A.cols());

    m_matrix->n = A.rows();
    m_matrix->nnz = A.nonZeros();

    m_matrix->fmttype = A.IsRowMajor ? SpmCSR : SpmCSC;
    m_matrix->layout = A.IsRowMajor ? SpmRowMajor : SpmColMajor;

    m_matrix->colptr = const_cast<pastix_int_t*>(A.IsRowMajor ? A.outerIndexPtr()
                                                              : A.innerIndexPtr());
    m_matrix->rowptr = const_cast<pastix_int_t*>(A.IsRowMajor ? A.innerIndexPtr()
                                                              : A.outerIndexPtr());
    m_matrix->values = const_cast<double*>(A.valuePtr());

    if (build_sparsity_pattern)
    {
        spmUpdateComputedFields(m_matrix.get());
        this->analyse_pattern();
        build_sparsity_pattern = false;
    }

    auto constexpr nrhs = 1;

    pastix_task_numfact(m_data, m_matrix.get());

    pastix_task_solve(m_data, nrhs, x.data(), m_matrix->n);

    x = b;

    pastix_task_refine(m_data,
                       m_matrix->n,
                       nrhs,
                       const_cast<double*>(b.data()),
                       m_matrix->n,
                       x.data(),
                       m_matrix->n);

    auto const end = std::chrono::steady_clock::now();
    std::chrono::duration<double> const elapsed_seconds = end - start;
    std::cout << std::string(6, ' ') << "PaStiX LU direct solver took " << elapsed_seconds.count()
              << "s\n";
}
}
