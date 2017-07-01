
#define CATCH_CONFIG_MAIN

#include <catch.hpp>

#include "mesh/BasicMesh.hpp"

#include "mesh/solid/MaterialCoordinates.hpp"
#include "mesh/solid/Submesh.hpp"
#include "mesh/solid/femMesh.hpp"

#include "assembler/solid/femMatrix.hpp"

#include <json/json.h>
#include <range/v3/view.hpp>

#include "CubeJson.hpp"

using namespace neon;

TEST_CASE("Nonlinear system solver test")
{
    using solid::MaterialCoordinates;
    using solid::femMesh;
    using solid::femMatrix;

    // Read in a cube mesh from the json input file and use this to
    // test the functionality of the basic mesh
    Json::Value mesh_data, material_data, simulation_data, solver_data;
    Json::Reader mesh_file, material_file, simulation_file, solver_file;

    REQUIRE(mesh_file.parse(cube_mesh_json().c_str(), mesh_data));
    REQUIRE(material_file.parse(material_data_json().c_str(), material_data));
    REQUIRE(simulation_file.parse(simulation_data_json().c_str(), simulation_data));
    REQUIRE(solver_file.parse(solver_data_json().c_str(), solver_data));

    // Create the test objects
    BasicMesh basic_mesh(mesh_data);
    NodalCoordinates nodal_coordinates(mesh_data);
    femMesh fem_mesh(basic_mesh, material_data, simulation_data);

    // Create the system
    femMatrix fem_matrix(fem_mesh, solver_data);

    fem_matrix.solve();

    // SECTION("Check initialisation") {}
}
