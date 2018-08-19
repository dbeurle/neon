
#include "conjugate_gradient_ocl.hpp"

#ifdef ENABLE_OCL

#include "exceptions.hpp"
#include "dmatrix_vector_product.hpp"
#include "numeric/float_compare.hpp"

#define VIENNACL_HAVE_EIGEN
#define VIENNACL_WITH_OPENCL

#include <viennacl/vector.hpp>
#include <viennacl/compressed_matrix.hpp>
#include <viennacl/linalg/jacobi_precond.hpp>
#include <viennacl/linalg/cg.hpp>

#include <cmath>
#include <iostream>

namespace neon
{
conjugate_gradient_ocl::conjugate_gradient_ocl() : iterative_linear_solver()
{
    viennacl::ocl::set_context_device_type(0, viennacl::ocl::gpu_tag());
}

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

conjugate_gradient_ocl::~conjugate_gradient_ocl() = default;

void conjugate_gradient_ocl::solve(sparse_matrix const& A, vector& x, vector const& b)
{
    using vcl_sparse_matrix = viennacl::compressed_matrix<double>;
    using vcl_vector = viennacl::vector<double>;

    viennacl::linalg::cg_tag solver_tag(residual_tolerance, max_iterations);

    vcl_sparse_matrix vcl_A(A.rows(), A.cols());
    vcl_vector vcl_b(b.rows());

    // Copy from Eigen objects to ViennaCL objects
    viennacl::copy(A, vcl_A);
    viennacl::copy(b, vcl_b);

    viennacl::linalg::jacobi_precond<vcl_sparse_matrix> vcl_jacobi(vcl_A,
                                                                   viennacl::linalg::jacobi_tag());

    // Conjugate gradient solver with preconditioner on GPU
    vcl_vector const vcl_x = viennacl::linalg::solve(vcl_A, vcl_b, solver_tag, vcl_jacobi);

    // Copy back to host
    viennacl::copy(vcl_x, x);

    viennacl::backend::finish();

    std::cout << std::string(6, ' ') << "Conjugate Gradient iterations: " << solver_tag.iters()
              << " (max. " << max_iterations << "), estimated error: " << solver_tag.error()
              << " (min. " << residual_tolerance << ")\n";
}
}
#endif
