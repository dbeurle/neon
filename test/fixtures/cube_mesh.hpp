
#pragma once

/// @file

#include <string>

// A cube mesh in the input mesh format
std::string json_cube_mesh();

std::string material_data_json();

std::string simulation_data_json();

std::string const nonlinear_options_disp_broken_json("{\"nonlinear_options\" : { "
                                                     "\"DisacementIncrementTolerance\" : 1.0e-5, "
                                                     "\"residual_tolerance\" : 0.001}}");

std::string const nonlinear_options_force_broken_json("{\"nonlinear_options\" : { "
                                                      "\"displacement_tolerance\" : "
                                                      "1.0e-5, "
                                                      "\"ResidForceTolerance\" : 0.001}}");

std::string simulation_data_traction_json();

std::string solver_data_json();
