
#pragma once

namespace neon
{
namespace mechanical
{
namespace solid
{
/** Alias for a solid mechanics type for a symmetric material model */
using constitutive_model = neon::constitutive_model<3, 6>;
}
namespace plane
{
/** Alias for a plane strain/stress type for a symmetric material model */
using constitutive_model = neon::constitutive_model<2, 3>;
}
}

namespace diffusion
{
/** Alias for a diffusion type for a general 3D material model */
using constitutive_model = neon::constitutive_model<3, 3>;
}
}
