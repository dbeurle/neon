
#include <catch2/catch.hpp>

#include "mesh/basic_mesh.hpp"
#include "mesh/material_coordinates.hpp"
#include "assembler/mechanics/latin_matrix.hpp"
#include "mesh/mechanics/solid/mesh.hpp"
#include "assembler/mechanics/static_matrix.hpp"
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
    using fem_mesh = neon::mechanics::solid::mesh<neon::mechanics::solid::submesh>;
    using static_matrix = neon::mechanics::static_matrix<fem_mesh>;

    // Read in a cube mesh from the json input file and use this to
    // test the functionality of the basic mesh
    neon::basic_mesh basic_mesh(json::parse(json_cube_mesh()));

    neon::nodal_coordinates nodal_coordinates(json::parse(json_cube_mesh()));

    auto simulation_data = json::parse(simulation_data_json());

    fem_mesh mesh(basic_mesh,
                  json::parse(material_data_json()),
                  simulation_data,
                  simulation_data["time"]["increments"]["initial"]);

    SECTION("Correct behaviour")
    {
        // Create the system and solve it
        static_matrix matrix(mesh, json::parse(simulation_data_json()));
        matrix.solve();
    }
}
TEST_CASE("LATIN solver test")
{
    using fem_mesh = neon::mechanics::solid::mesh<neon::mechanics::solid::latin_submesh>;
    using latin_matrix = neon::mechanics::latin_matrix<fem_mesh>;

    // Read in a cube mesh from the json input file and use this to
    // test the functionality of the basic mesh
    neon::basic_mesh basic_mesh(json::parse(json_cube_mesh()));

    neon::nodal_coordinates nodal_coordinates(json::parse(json_cube_mesh()));

    auto simulation_data = json::parse(simulation_data_json());

    fem_mesh mesh(basic_mesh,
                  json::parse(material_data_json()),
                  simulation_data,
                  simulation_data["time"]["increments"]["initial"]);

    SECTION("Correct behaviour")
    {
        // Create the system and solve it
        latin_matrix matrix(mesh, json::parse(simulation_data_json()));
        matrix.solve();
    }
}
