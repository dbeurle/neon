
#include "eigen_solver.hpp"
#include "Exceptions.hpp"

#include <iostream>
#include "arpack_ng/ArpackSelfAdjointEigenSolver.h"

namespace neon
{
eigen_solver::eigen_solver(std::int64_t const eigenvalues_to_extract)
    : eigenvalues_to_extract{eigenvalues_to_extract}
{
}

std::pair<vector, matrix> eigen_solver::solve(sparse_matrix const& A)
{
    Eigen::SparseMatrix<double> A_col = A;

    Eigen::ArpackGeneralizedSelfAdjointEigenSolver<decltype(A_col)> arpack;

    arpack.compute(A_col, eigenvalues_to_extract, "LM");

    if (arpack.getNbrConvergedEigenValues() < eigenvalues_to_extract)
    {
        throw computational_error("Eigenvalues did not converge");
    }

    if (arpack.info() != Eigen::Success)
    {
        throw computational_error("Numerical issued occurred");
    }

    return {arpack.eigenvectors(), arpack.eigenvalues()};
}

std::pair<vector, matrix> eigen_solver::solve(sparse_matrix const& A, sparse_matrix const& B)
{
    Eigen::SparseMatrix<double> A_col = A;
    Eigen::SparseMatrix<double> B_col = B;

    Eigen::ArpackGeneralizedSelfAdjointEigenSolver<decltype(A_col)> arpack;

    arpack.compute(A_col, B_col, eigenvalues_to_extract, "LM");

    if (arpack.getNbrConvergedEigenValues() < eigenvalues_to_extract)
    {
        throw computational_error("Eigenvalues did not converge");
    }

    if (arpack.info() != Eigen::Success)
    {
        throw computational_error("Numerical issued occurred");
    }

    return {arpack.eigenvectors(), arpack.eigenvalues()};
}
}
