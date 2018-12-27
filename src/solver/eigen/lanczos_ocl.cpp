
#include "solver/eigen/lanczos_ocl.hpp"

#ifdef ENABLE_OPENCL

#define VIENNACL_HAVE_EIGEN

#include <viennacl/compressed_matrix.hpp>
#include <viennacl/matrix.hpp>
#include <viennacl/matrix_proxy.hpp>
#include <viennacl/vector.hpp>

#include <viennacl/linalg/lanczos.hpp>

namespace neon
{
lanczos_ocl::lanczos_ocl(std::int64_t const values_to_extract,
                         eigen_solver::eigen_spectrum const spectrum)
    : eigen_solver(values_to_extract, spectrum)
{
}

void lanczos_ocl::solve(sparse_matrix const& A)
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

void lanczos_ocl::solve(sparse_matrix const&, sparse_matrix const&) {}
}
#endif
