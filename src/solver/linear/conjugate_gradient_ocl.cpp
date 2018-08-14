
#include "conjugate_gradient_ocl.hpp"

#ifdef ENABLE_OCL

#include "exceptions.hpp"
#include "dmatrix_vector_product.hpp"
#include "numeric/float_compare.hpp"

#define VIENNACL_HAVE_EIGEN 1

#include <viennacl/vector.hpp>
#include <viennacl/compressed_matrix.hpp>
#include <viennacl/linalg/ilu.hpp>
#include <viennacl/linalg/jacobi_precond.hpp>
#include <viennacl/linalg/cg.hpp>

#include <cmath>
#include <iostream>

namespace neon
{
conjugate_gradient_ocl::conjugate_gradient_ocl() : iterative_linear_solver() {}

conjugate_gradient_ocl::conjugate_gradient_ocl(double const residual_tolerance)
    : conjugate_gradient_ocl()
{
    this->residual_tolerance = residual_tolerance;
}

conjugate_gradient_ocl::conjugate_gradient_ocl(int const max_iterations) : conjugate_gradient_ocl()
{
    this->max_iterations = max_iterations;
}

conjugate_gradient_ocl::conjugate_gradient_ocl(double const residual_tolerance,
                                               int const max_iterations)
    : conjugate_gradient_ocl()
{
    this->residual_tolerance = residual_tolerance;
    this->max_iterations = max_iterations;
}

conjugate_gradient_ocl::~conjugate_gradient_ocl() {}

void conjugate_gradient_ocl::solve(sparse_matrix const& A, vector& x, vector const& b)
{
    using vcl_sparse_matrix = viennacl::compressed_matrix<double>;
    using vcl_vector = viennacl::vector<double>;

    vcl_sparse_matrix vcl_A(A.rows(), A.cols());
    vcl_vector vcl_x(x.size());
    vcl_vector vcl_b(b.size());

    // Copy from Eigen objects to ViennaCL objects
    viennacl::copy(A, vcl_A);
    viennacl::copy(x, vcl_x);
    viennacl::copy(b, vcl_b);

    // Compute ILUT preconditioners for CPU and for GPU objects:
    // 10 entries, rel. tol. 1e-5
    viennacl::linalg::ilut_tag ilut_conf(10, 1e-5);

    // preconditioner for ViennaCL objects:
    viennacl::linalg::ilut_precond<vcl_sparse_matrix> vcl_ilut(vcl_A, ilut_conf);

    // Conjugate gradient solver without preconditioner on GPU
    vcl_x = viennacl::linalg::solve(vcl_A, vcl_b, viennacl::linalg::cg_tag());
    // Conjugate gradient solver using ILUT preconditioner on GPU
    vcl_x = viennacl::linalg::solve(vcl_A, vcl_b, viennacl::linalg::cg_tag(), vcl_ilut);

    // Copy back to host
    viennacl::copy(vcl_x, x);
}
}
#endif
