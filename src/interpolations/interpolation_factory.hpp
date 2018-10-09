
#pragma once

#include "shape_function_forward.hpp"
#include "mesh/element_topology.hpp"
#include "math/bilinear_form.hpp"
#include "io/json_forward.hpp"

#include <memory>

/// \file interpolation_factory.hpp

namespace neon
{
/// Factory method for the three dimensional shape functions
/// \param topology Element topology
/// \param simulation_data json object for element technology input options
std::unique_ptr<volume_interpolation> make_volume_interpolation(
    element_topology const topology,
    json const& simulation_data,
    std::variant<fem::linear, fem::bilinear> integral_form,
    std::optional<fem::bilinear> bilinear_form);

/// Factory method for the two dimensional shape functions
/// \param topology Element topology
/// \param simulation_data json object for element technology input options
std::unique_ptr<surface_interpolation> make_surface_interpolation(element_topology const topology,
                                                                  json const& simulation_data);

/// Factory method for the two dimensional shape functions
/// \param topology Element topology
/// \param simulation_data json object for element technology input options
std::unique_ptr<line_interpolation> make_line_interpolation(element_topology const topology,
                                                            json const& simulation_data);
}
