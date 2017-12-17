
#pragma once

#include "assembler/mechanical/femStaticMatrix.hpp"

#include "mesh/mechanical/plane/femMesh.hpp"

namespace neon::mechanical::plane
{
using femStaticMatrix = detail::femStaticMatrix<plane::femMesh>;
}
