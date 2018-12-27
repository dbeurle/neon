
#pragma once

#include "shape_function_forward.hpp"
#include "mesh/element_topology.hpp"

#include <memory>

/// \file interpolation_factory.hpp

namespace neon
{
/// Factory method for the three dimensional shape functions
/// \param topology Element topology
std::unique_ptr<volume_interpolation> make_volume_interpolation(element_topology const topology);

/// Factory method for the two dimensional shape functions
/// \param topology Element topology
std::unique_ptr<surface_interpolation> make_surface_interpolation(element_topology const topology);

/// Factory method for the two dimensional shape functions
/// \param topology Element topology
std::unique_ptr<line_interpolation> make_line_interpolation(element_topology const topology);
}
