
#pragma once

#include "constitutive/VoigtDimension.hpp"

namespace neon
{
template <int rank2_dimension, int rank4_dimension>
class InternalVariables;

namespace mechanical
{
namespace solid
{
using InternalVariables = neon::InternalVariables<3, 6>;
}
namespace plane
{
using InternalVariables = neon::InternalVariables<2, 3>;
}
}

namespace diffusion
{
using InternalVariables = neon::InternalVariables<3, 3>;
}
}
