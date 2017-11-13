
#define CATCH_CONFIG_MAIN

#include <catch.hpp>

#include "mesh/BasicMesh.hpp"

#include "mesh/MaterialCoordinates.hpp"
#include "mesh/solid/femSubmesh.hpp"
#include "mesh/solid/femMesh.hpp"

#include "assembler/solid/femStaticMatrix.hpp"

#include "assembler/solid/femDynamicMatrix.hpp"
#include "solver/TimeControl.hpp"

#include <json/json.h>
#include <range/v3/view.hpp>

#include "CubeJson.hpp"

using namespace neon;

TEST_CASE("Nonlinear system equilibrium solver test")
{
    using mech::solid::femMesh;
    using mech::solid::femStaticMatrix;

    // Read in a cube mesh from the json input file and use this to
    // test the functionality of the basic mesh
    Json::Value mesh_data, material_data, simulation_data;
    Json::Reader reader;

    REQUIRE(reader.parse(cube_mesh_json().c_str(), mesh_data));
    REQUIRE(reader.parse(material_data_json().c_str(), material_data));
    REQUIRE(reader.parse(simulation_data_json().c_str(), simulation_data));

    // Create the test objects
    BasicMesh basic_mesh(mesh_data);
    NodalCoordinates nodal_coordinates(mesh_data);
    femMesh fem_mesh(basic_mesh, material_data, simulation_data);

    SECTION("Correct behaviour")
    {
        // Create the system
        femStaticMatrix fem_matrix(fem_mesh, simulation_data);

        fem_matrix.solve();
    }
    // SECTION("NonlinearOptions displacement increment incorrect")
    // {
    //     REQUIRE(reader.parse(nonlinear_options_disp_broken_json.c_str(), nonlinear_data));
    //
    //     REQUIRE_THROWS_AS(femStaticMatrix(fem_mesh, simulation_data), std::runtime_error);
    // }
    // SECTION("NonlinearOptions residual force incorrect")
    // {
    //     REQUIRE(reader.parse(nonlinear_options_force_broken_json.c_str(), nonlinear_data));
    //
    //     REQUIRE_THROWS_AS(femStaticMatrix(fem_mesh, sosimulation_data), std::runtime_error);
    // }
}
