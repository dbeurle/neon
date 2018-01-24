
#define CATCH_CONFIG_MAIN

#include <catch.hpp>

#include "mesh/BasicMesh.hpp"

#include "mesh/MaterialCoordinates.hpp"
#include "mesh/mechanical/solid/femMesh.hpp"
#include "mesh/mechanical/solid/femSubmesh.hpp"

#include "assembler/mechanical/solid/femStaticMatrix.hpp"

#include "solver/TimeControl.hpp"

#include "io/json.hpp"

#include <range/v3/view.hpp>

#include "CubeJson.hpp"

using namespace neon;

TEST_CASE("Nonlinear system equilibrium solver test")
{
    using mechanical::solid::femMesh;
    using mechanical::solid::femStaticMatrix;

    // Read in a cube mesh from the json input file and use this to
    // test the functionality of the basic mesh
    BasicMesh basic_mesh(json::parse(cube_mesh_json()));

    NodalCoordinates nodal_coordinates(json::parse(cube_mesh_json()));

    femMesh fem_mesh(basic_mesh,
                     json::parse(material_data_json()),
                     json::parse(simulation_data_json()));

    SECTION("Correct behaviour")
    {
        // Create the system
        femStaticMatrix fem_matrix(fem_mesh, json::parse(simulation_data_json()));

        fem_matrix.solve();
    }
}
