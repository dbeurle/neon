
#pragma once

#include "constitutive/VoigtDimension.hpp"

namespace neon
{
template <int spatial_dimension, int voigt_dimension = spatial_to_voigt(spatial_dimension)>
class InternalVariables;

namespace mechanical
{
namespace solid
{
using InternalVariables = neon::InternalVariables<3>;
}
namespace plane
{
using InternalVariables = neon::InternalVariables<2>;
}
}

namespace diffusion
{
using InternalVariables = neon::InternalVariables<3>;
}
}
