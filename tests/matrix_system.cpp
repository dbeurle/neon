
#define CATCH_CONFIG_MAIN

#include <catch.hpp>

#include "mesh/basic_mesh.hpp"

#include "mesh/material_coordinates.hpp"
#include "mesh/mechanical/solid/fem_mesh.hpp"
#include "mesh/mechanical/solid/fem_submesh.hpp"

#include "assembler/mechanical/solid/fem_static_matrix.hpp"

#include "solver/time_step_control.hpp"

#include "io/json.hpp"

#include <range/v3/view.hpp>

#include "CubeJson.hpp"

using namespace neon;

TEST_CASE("Nonlinear system equilibrium solver test")
{
    using mechanical::solid::fem_mesh;
    using mechanical::solid::fem_static_matrix;

    // Read in a cube mesh from the json input file and use this to
    // test the functionality of the basic mesh
    basic_mesh basic_mesh(json::parse(cube_mesh_json()));

    nodal_coordinates nodal_coordinates(json::parse(cube_mesh_json()));

    fem_mesh fem_mesh(basic_mesh,
                     json::parse(material_data_json()),
                     json::parse(simulation_data_json()));

    SECTION("Correct behaviour")
    {
        // Create the system
        fem_static_matrix fem_matrix(fem_mesh, json::parse(simulation_data_json()));

        fem_matrix.solve();
    }
}
