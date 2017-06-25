/*
 * neon - A finite element solver.
 *
 * For licensing please refer to the LICENSE.md file
 *
 */

#include "MaterialCoordinates.hpp"

namespace neon::solid
{
MaterialCoordinates::MaterialCoordinates(Vector const& initial_coordinates)
    : NodalCoordinates(initial_coordinates), x(initial_coordinates)
{
}
}
