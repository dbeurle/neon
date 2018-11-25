
#include "eigenvalue_solver.hpp"
#include "exceptions.hpp"

#include "arpack_ng/ArpackSelfAdjointEigenSolver.h"

#define VIENNACL_HAVE_EIGEN

#include <viennacl/compressed_matrix.hpp>
#include <viennacl/matrix.hpp>
#include <viennacl/matrix_proxy.hpp>
#include <viennacl/vector.hpp>

#include <viennacl/linalg/power_iter.hpp>
#include <viennacl/linalg/lanczos.hpp>

#include <vector>

namespace neon
{
eigenvalue_solver::eigenvalue_solver(std::int64_t const values_to_extract)
    : values_to_extract{values_to_extract}
{
}

void eigenvalue_solver::solve(sparse_matrix const& A)
{
    Eigen::SparseMatrix<double> A_col = A;

    Eigen::ArpackGeneralizedSelfAdjointEigenSolver<decltype(A_col)> arpack;

    arpack.compute(A_col, values_to_extract, "SM");

    if (arpack.getNbrConvergedEigenValues() < values_to_extract)
    {
        throw computational_error("Eigenvalues did not converge");
    }

    if (arpack.info() != Eigen::Success)
    {
        throw computational_error("Numerical issued occurred");
    }

    m_eigenvalues = arpack.eigenvectors();
    m_eigenvectors = arpack.eigenvalues();
}

void eigenvalue_solver::solve(sparse_matrix const& A, sparse_matrix const& B)
{
    Eigen::SparseMatrix<double> A_col = A;
    Eigen::SparseMatrix<double> B_col = B;

    Eigen::ArpackGeneralizedSelfAdjointEigenSolver<decltype(A_col)> arpack;

    arpack.compute(A_col, B_col, values_to_extract, "SM");

    if (arpack.getNbrConvergedEigenValues() < values_to_extract)
    {
        throw computational_error("Eigenvalues did not converge");
    }

    if (arpack.info() != Eigen::Success)
    {
        throw computational_error("Numerical issued occurred");
    }

    m_eigenvalues = arpack.eigenvalues();
    m_eigenvectors = arpack.eigenvectors();
}

power_iteration::power_iteration(std::int64_t const values_to_extract) : eigenvalue_solver(1) {}

void power_iteration::solve(sparse_matrix const& A)
{
    m_eigenvalues.resize(1);
    m_eigenvectors.resize(A.rows(), 1);

    viennacl::linalg::power_iter_tag power_tag(1.0e-8);

    viennacl::compressed_matrix<double> vcl_A(A.rows(), A.cols());
    viennacl::vector<double> vcl_eigenvectors(A.rows());

    viennacl::copy(A, vcl_A);

    m_eigenvalues(0) = viennacl::linalg::eig(vcl_A, power_tag, vcl_eigenvectors);

    // We can't copy into the first column with Eigen, so take a copy
    vector eigenvector(A.rows());
    viennacl::copy(vcl_eigenvectors, eigenvector);

    m_eigenvectors.col(0) = eigenvector.normalized();
}

void power_iteration::solve(sparse_matrix const& A, sparse_matrix const& B)
{
    throw std::runtime_error("Method not implemented " + std::string(__FUNCTION__));
}

lanzcos_solver::lanzcos_solver(std::int64_t const values_to_extract)
    : eigenvalue_solver(values_to_extract)
{
}

void lanzcos_solver::solve(sparse_matrix const& A)
{
    viennacl::linalg::lanczos_tag lanczos_tag(0.85,
                                              values_to_extract,
                                              viennacl::linalg::lanczos_tag::partial_reorthogonalization,
                                              200);

    viennacl::compressed_matrix<double> vcl_A(A.rows(), A.cols());
    viennacl::matrix<double> vcl_eigenvectors(A.rows(), lanczos_tag.num_eigenvalues());

    m_eigenvalues.resize(lanczos_tag.num_eigenvalues());
    m_eigenvectors.resize(A.rows(), lanczos_tag.num_eigenvalues());

    viennacl::copy(A, vcl_A);

    auto const vcl_eigenvalues = viennacl::linalg::eig(vcl_A, vcl_eigenvectors, lanczos_tag);

    std::copy(begin(vcl_eigenvalues), end(vcl_eigenvalues), m_eigenvalues.data());
    viennacl::copy(vcl_eigenvectors, m_eigenvectors);

    for (std::int64_t column{}; column < values_to_extract; ++column)
    {
        m_eigenvectors.col(column) = m_eigenvectors.col(column).normalized().eval();
    }
}

void lanzcos_solver::solve(sparse_matrix const& A, sparse_matrix const& B) {}
}
