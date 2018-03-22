
#pragma once

#include "numeric/dense_matrix.hpp"

/// \file Projection.hpp

namespace neon::geometry
{
matrix2x project_to_plane(matrix3x const& nodal_coordinates);

vector3 unit_outward_normal(matrix3x const& nodal_coordinates);
}
