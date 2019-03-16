
#pragma once

#include "numeric/dense_matrix.hpp"

/// @file

namespace neon
{
namespace mechanics
{
namespace solid
{
/// Alias for a solid mechanics type for a symmetric material model
using constitutive_model = neon::constitutive_model<matrix3, matrix6>;
}
namespace plane
{
/// Alias for a plane strain/stress type for a symmetric material model
using constitutive_model = neon::constitutive_model<matrix2, matrix3>;
}
}

namespace diffusion
{
/// Alias for a diffusion type for a general 3D material model
using constitutive_model = neon::constitutive_model<matrix3, matrix3>;

namespace reaction
{
using constitutive_model = neon::constitutive_model<matrix3, matrix3>;
}
}
}
