
#include <catch.hpp>

#include "mesh/basic_mesh.hpp"
#include "mesh/material_coordinates.hpp"
#include "mesh/mechanical/solid/fem_mesh.hpp"
#include "assembler/mechanical/fem_static_matrix.hpp"
#include "io/json.hpp"

#include "CubeJson.hpp"

TEST_CASE("Nonlinear system equilibrium solver test")
{
    using fem_mesh = neon::mechanical::solid::fem_mesh;
    using fem_static_matrix = neon::mechanical::fem_static_matrix<fem_mesh>;

    // Read in a cube mesh from the json input file and use this to
    // test the functionality of the basic mesh
    neon::basic_mesh basic_mesh(json::parse(cube_mesh_json()));

    neon::nodal_coordinates nodal_coordinates(json::parse(cube_mesh_json()));

    auto simulation_data = json::parse(simulation_data_json());

    fem_mesh mesh(basic_mesh,
                  json::parse(material_data_json()),
                  simulation_data,
                  simulation_data["Time"]["Increments"]["Initial"]);

    SECTION("Correct behaviour")
    {
        // Create the system and solve it
        fem_static_matrix matrix(mesh, json::parse(simulation_data_json()));
        matrix.solve();
    }
}
