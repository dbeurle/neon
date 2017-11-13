
#pragma once

#include "constitutive/VoigtDimension.hpp"

namespace neon
{
template <int spatial_dimension,
          int voigt_dimension = spatial_to_voigt(std::integral_constant<int, spatial_dimension>{})>
class InternalVariables;

namespace mech::solid
{
using InternalVariables = neon::InternalVariables<3>;
}
namespace mech::beam
{
using InternalVariables = neon::InternalVariables<1>;
}

namespace diffusion
{
using InternalVariables = neon::InternalVariables<3>;
}
}
