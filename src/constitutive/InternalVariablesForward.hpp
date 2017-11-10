
#pragma once

namespace neon
{
template <int spatial_dimension>
class InternalVariables;

namespace solid
{
using InternalVariables = neon::InternalVariables<3>;
}
namespace beam
{
using InternalVariables = neon::InternalVariables<1>;
}

namespace diffusion
{
using InternalVariables = neon::InternalVariables<3>;
}
}
