
#include "solver/eigen/arpack.hpp"
#include "exceptions.hpp"

#include "arpack_ng/ArpackSelfAdjointEigenSolver.h"

namespace neon
{
arpack::arpack(std::int64_t const values_to_extract, eigen_solver::eigen_spectrum const spectrum)
    : eigen_solver{values_to_extract, spectrum}
{
}

void arpack::solve(sparse_matrix const& A)
{
    Eigen::SparseMatrix<double> A_col = A;

    Eigen::ArpackGeneralizedSelfAdjointEigenSolver<decltype(A_col)> arpack;

    arpack.compute(A_col,
                   values_to_extract,
                   (m_spectrum == eigen_solver::eigen_spectrum::lower ? "SM" : "LM"));

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

void arpack::solve(sparse_matrix const& A, sparse_matrix const& B)
{
    Eigen::SparseMatrix<double> A_col = A;
    Eigen::SparseMatrix<double> B_col = B;

    Eigen::ArpackGeneralizedSelfAdjointEigenSolver<decltype(A_col)> arpack;

    arpack.compute(A_col,
                   B_col,
                   values_to_extract,
                   (m_spectrum == eigen_solver::eigen_spectrum::lower ? "SM" : "LM"));

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
}
