
#pragma once

/// @file

namespace neon
{
namespace mechanics
{
namespace solid
{
using internal_variables_t = neon::internal_variables<3, 6>;
}
namespace plane
{
using internal_variables_t = neon::internal_variables<2, 3>;
}
namespace beam
{
using internal_variables_t = neon::internal_variables<2, 2>;
}
}

namespace diffusion
{
using internal_variables_t = neon::internal_variables<3, 3>;
}
}
