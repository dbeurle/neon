
#include "solver/linear/biconjugate_gradient_stabilised_ocl.hpp"

#ifdef ENABLE_OCL

#include "exceptions.hpp"
#include "dmatrix_vector_product.hpp"
#include "numeric/float_compare.hpp"

#include <cmath>
#include <iostream>

namespace neon
{
void biconjugate_gradient_stabilised_ocl::solve(sparse_matrix const& A, vector& x, vector const& b)
{
}
}
#endif
