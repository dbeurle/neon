
#pragma once

/// @file

#include <cstdint>

namespace neon
{
namespace cuda
{
void diagonal_matrix_vector_product(double const* const diagonal_matrix,
                                    double* const vector,
                                    std::size_t const size);
}
}
