
#pragma once

#include "assembler/mechanical/fem_static_matrix.hpp"
#include "mesh/mechanical/solid/fem_mesh.hpp"

namespace neon::mechanical::solid
{
using fem_static_matrix = detail::fem_static_matrix<solid::fem_mesh>;
}
