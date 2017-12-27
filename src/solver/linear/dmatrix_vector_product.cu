
#include "dmatrix_vector_product.hpp"

namespace neon
{
namespace cuda
{
__global__ void diagonal_matrix_vector_product_kernel(double const* const __restrict__ diagonal_matrix,
                                                      double* const __restrict__ vector,
                                                      std::size_t const size)
{
    auto const i = blockIdx.x * blockDim.x + threadIdx.x;

    if (i >= size) return;

    vector[i] *= diagonal_matrix[i];
}

void diagonal_matrix_vector_product(double const* const diagonal_matrix,
                                    double* const vector,
                                    std::size_t const size)
{
    diagonal_matrix_vector_product_kernel<<<size / 256 + 1, 256>>>(diagonal_matrix, vector, size);
}
}
}
