
#pragma once

#include "numeric/dense_matrix.hpp"

/// @file

namespace neon
{
namespace mechanics
{
namespace solid
{
using internal_variables_t = neon::internal_variables<matrix3, matrix6>;
}
namespace plane
{
using internal_variables_t = neon::internal_variables<matrix2, matrix3>;
}
namespace beam
{
using internal_variables_t = neon::internal_variables<matrix2, matrix2>;
}
}

namespace diffusion
{
using internal_variables_t = neon::internal_variables<matrix3, matrix3>;
}
}
