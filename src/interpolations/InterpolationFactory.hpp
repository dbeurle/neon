
#pragma once

#include "ShapeFunction.hpp"

#include "mesh/ElementTopology.hpp"

#include <json/forwards.h>
#include <memory>

namespace neon
{
namespace solid
{
/** Factory method for the three dimensional shape functions */
std::unique_ptr<VolumeInterpolation> make_shape_function(ElementTopology const topology,
                                                         Json::Value const& simulation_data);
}
}
