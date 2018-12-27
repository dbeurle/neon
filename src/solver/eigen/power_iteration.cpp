
#include "solver/eigen/power_iteration.hpp"

#ifdef ENABLE_OPENCL

#include "exceptions.hpp"

#define VIENNACL_HAVE_EIGEN

#include <viennacl/compressed_matrix.hpp>
#include <viennacl/matrix.hpp>
#include <viennacl/matrix_proxy.hpp>
#include <viennacl/vector.hpp>

#include <viennacl/linalg/power_iter.hpp>

namespace neon
{
power_iteration::power_iteration(std::int64_t, eigen_solver::eigen_spectrum const spectrum)
    : eigen_solver(1, spectrum)
{
}

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

void power_iteration::solve(sparse_matrix const&, sparse_matrix const&)
{
    throw std::runtime_error("Method not implemented " + std::string(__FUNCTION__));
}
}
#endif
