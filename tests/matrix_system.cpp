
#include <catch.hpp>

#include "mesh/basic_mesh.hpp"
#include "mesh/material_coordinates.hpp"
#include "mesh/mechanical/solid/fem_mesh.hpp"
#include "assembler/mechanical/fem_static_matrix.hpp"
#include "numeric/doublet.hpp"
#include "io/json.hpp"

#include "fixtures/cube_mesh.hpp"

using neon::json;

TEST_CASE("Doublet class")
{
    SECTION("Trival zero case")
    {
        neon::doublet row_col_val(0, 0);

        REQUIRE(row_col_val.row() == 0);
        REQUIRE(row_col_val.col() == 0);
        REQUIRE(row_col_val.value() == Approx(0.0));
    }
    SECTION("Nontrivial case")
    {
        neon::doublet row_col_val(3, 1);

        REQUIRE(row_col_val.row() == 3);
        REQUIRE(row_col_val.col() == 1);
        REQUIRE(row_col_val.value() == Approx(0.0));
    }
    SECTION("Nontrivial case reversed")
    {
        neon::doublet row_col_val(1, 3);

        REQUIRE(row_col_val.row() == 1);
        REQUIRE(row_col_val.col() == 3);
        REQUIRE(row_col_val.value() == Approx(0.0));
    }
}
TEST_CASE("Nonlinear system equilibrium solver test")
{
    using fem_mesh = neon::mechanical::solid::fem_mesh;
    using fem_static_matrix = neon::mechanical::fem_static_matrix<fem_mesh>;

    // Read in a cube mesh from the json input file and use this to
    // test the functionality of the basic mesh
    neon::basic_mesh basic_mesh(json::parse(json_cube_mesh()));

    neon::nodal_coordinates nodal_coordinates(json::parse(json_cube_mesh()));

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
