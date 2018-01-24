
#pragma once

#include "shape_function.hpp"

#include "mesh/element_topology.hpp"

// TODO forward declaration
#include "io/json.hpp"

#include <memory>

namespace neon
{
/** Factory method for the three dimensional shape functions */
std::unique_ptr<volume_interpolation> make_volume_interpolation(element_topology const topology,
                                                                json const& simulation_data);

/** Factory method for the two dimensional shape functions */
std::unique_ptr<surface_interpolation> make_surface_interpolation(element_topology const topology,
                                                                  json const& simulation_data);

/** Factory method for the two dimensional shape functions */
std::unique_ptr<line_interpolation> make_line_interpolation(element_topology const topology,
                                                            json const& simulation_data);
}
