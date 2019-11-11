
#include "PaStiX.hpp"

#include <chrono>
#include <termcolor/termcolor.hpp>

#include <unsupported/Eigen/SparseExtra>

namespace neon
{
PaStiX::PaStiX() : m_matrix(std::make_unique<spmatrix_t>())
{
    pastixInitParam(m_integer_parameters.data(), m_float_parameters.data());

    m_integer_parameters[IPARM_FACTORIZATION] = PastixFactLDLT;
    m_integer_parameters[IPARM_VERBOSE] = PastixVerboseNo;
    m_integer_parameters[IPARM_THREAD_NBR] = 1;

    pastixInit(&m_data, MPI_COMM_WORLD, m_integer_parameters.data(), m_float_parameters.data());

    m_matrix->mtxtype = SpmGeneral;
    m_matrix->flttype = std::is_same<sparse_matrix::value_type, double>::value ? SpmDouble : SpmFloat;
    m_matrix->fmttype = SpmCSR;

    m_matrix->gN = 0;
    m_matrix->n = 0;
    m_matrix->gnnz = 0;
    m_matrix->nnz = 0;

    m_matrix->gNexp = 0;
    m_matrix->nexp = 0;
    m_matrix->gnnzexp = 0;
    m_matrix->nnzexp = 0;

    m_matrix->dof = 1;
    m_matrix->dofs = nullptr;
    m_matrix->layout = SpmRowMajor;

    m_matrix->colptr = nullptr;
    m_matrix->rowptr = nullptr;
    m_matrix->loc2glob = nullptr;
    m_matrix->values = nullptr;
}

PaStiX::~PaStiX() { pastixFinalize(&m_data); }

void PaStiX::analyse_pattern()
{
    spmUpdateComputedFields(m_matrix.get());

    if (auto const code = pastix_task_analyze(m_data, m_matrix.get()); code != PASTIX_SUCCESS)
    {
        throw std::runtime_error("PaStiX: Analyze routine exited with code " + std::to_string(code));
    }
}

PaStiXLDLT::PaStiXLDLT() { m_integer_parameters[IPARM_FACTORIZATION] = PastixFactLLT; }

void PaStiXLDLT::solve(sparse_matrix const& A, vector& x, vector const& b)
{
    auto const start = std::chrono::steady_clock::now();

    m_matrix->colptr = const_cast<pastix_int_t*>(A.IsRowMajor ? A.outerIndexPtr()
                                                              : A.innerIndexPtr());
    m_matrix->rowptr = const_cast<pastix_int_t*>(A.IsRowMajor ? A.innerIndexPtr()
                                                              : A.outerIndexPtr());
    m_matrix->values = const_cast<double*>(A.valuePtr());

    x = b;

    if (build_sparsity_pattern)
    {
        m_matrix->n = A.rows();
        m_matrix->nnz = A.nonZeros();

        m_matrix->fmttype = A.IsRowMajor ? SpmCSR : SpmCSC;
        m_matrix->layout = A.IsRowMajor ? SpmRowMajor : SpmColMajor;

        this->analyse_pattern();
        build_sparsity_pattern = false;
    }

    auto constexpr nrhs = 1;

    if (auto const code = pastix_task_numfact(m_data, m_matrix.get()); code != PASTIX_SUCCESS)
    {
        throw std::runtime_error("PaStiXLDLT: Factorisation routine exited with code "
                                 + std::to_string(code));
    }

    if (auto const code = pastix_task_solve(m_data, nrhs, x.data(), m_matrix->n);
        code != PASTIX_SUCCESS)
    {
        throw std::runtime_error("PaStiXLDLT: Solve routine exited with code " + std::to_string(code));
    }

    if (auto const code = pastix_task_refine(m_data,
                                             m_matrix->n,
                                             nrhs,
                                             const_cast<double*>(b.data()),
                                             m_matrix->n,
                                             x.data(),
                                             m_matrix->n);
        code != PASTIX_SUCCESS)
    {
        throw std::runtime_error("PaStiXLDLT: Refine routine exited with code "
                                 + std::to_string(code));
    }

    auto const end = std::chrono::steady_clock::now();
    std::chrono::duration<double> const elapsed_seconds = end - start;
    std::cout << std::string(6, ' ') << "PaStiX LDLT direct solver took " << elapsed_seconds.count()
              << "s\n";
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
