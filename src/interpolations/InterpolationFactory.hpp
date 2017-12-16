
#pragma once

#include "ShapeFunction.hpp"

#include "mesh/ElementTopology.hpp"

#include <json/forwards.h>
#include <memory>

namespace neon
{
/** Factory method for the three dimensional shape functions */
std::unique_ptr<VolumeInterpolation> make_volume_interpolation(ElementTopology const topology,
                                                               Json::Value const& simulation_data);

/** Factory method for the two dimensional shape functions */
std::unique_ptr<SurfaceInterpolation> make_surface_interpolation(ElementTopology const topology,
                                                                 Json::Value const& simulation_data);

/** Factory method for the two dimensional shape functions */
std::unique_ptr<LineInterpolation> make_line_interpolation(ElementTopology const topology,
                                                           Json::Value const& simulation_data);
}
