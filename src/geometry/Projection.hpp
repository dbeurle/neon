
#pragma once

#include "numeric/DenseMatrix.hpp"

namespace neon::geometry
{
Matrix2x project_to_plane(Matrix3x const& nodal_coordinates);

Vector3 unit_outward_normal(Matrix3x const& nodal_coordinates);
}
