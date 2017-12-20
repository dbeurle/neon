
#define CATCH_CONFIG_MAIN

#include <catch.hpp>

#include "mesh/BasicMesh.hpp"

#include "mesh/MaterialCoordinates.hpp"
#include "mesh/mechanical/solid/femMesh.hpp"
#include "mesh/mechanical/solid/femSubmesh.hpp"

#include "assembler/solid/femStaticMatrix.hpp"

#include "assembler/solid/femDynamicMatrix.hpp"
#include "solver/TimeControl.hpp"

#include <json/json.h>
#include <range/v3/view.hpp>

#include "CubeJson.hpp"

Json::CharReaderBuilder reader;
JSONCPP_STRING input_errors;

using namespace neon;

TEST_CASE("Nonlinear system equilibrium solver test")
{
    using mechanical::solid::femMesh;
    using mechanical::solid::femStaticMatrix;

    // Read in a cube mesh from the json input file and use this to
    // test the functionality of the basic mesh
    Json::Value mesh_data, material_data, simulation_data;

    std::istringstream cube_mesh_stream(cube_mesh_json());
    std::istringstream material_data_stream(material_data_json());
    std::istringstream simulation_data_stream(simulation_data_json());

    REQUIRE(Json::parseFromStream(reader, cube_mesh_stream, &mesh_data, &input_errors));
    REQUIRE(Json::parseFromStream(reader, material_data_stream, &material_data, &input_errors));
    REQUIRE(Json::parseFromStream(reader, simulation_data_stream, &simulation_data, &input_errors));

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
}
