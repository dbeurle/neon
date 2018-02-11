
#pragma once

#include "assembler/mechanical/femStaticMatrix.hpp"

#include "mesh/mechanical/plane/fem_mesh.hpp"

namespace neon::mechanical::plane
{
using femStaticMatrix = detail::femStaticMatrix<plane::fem_mesh>;
}
