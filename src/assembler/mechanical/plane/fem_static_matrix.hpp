
#pragma once

#include "assembler/mechanical/fem_static_matrix.hpp"
#include "mesh/mechanical/plane/fem_mesh.hpp"

namespace neon::mechanical::plane
{
using fem_static_matrix = detail::fem_static_matrix<plane::fem_mesh>;
}
