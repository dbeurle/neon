
#define CATCH_CONFIG_MAIN

#include <catch.hpp>

#include "mesh/BasicMesh.hpp"

#include "mesh/solid/MaterialCoordinates.hpp"
#include "mesh/solid/Submesh.hpp"
#include "mesh/solid/femMesh.hpp"

#include "assembler/solid/femStaticMatrix.hpp"

#include "assembler/solid/femDynamicMatrix.hpp"
#include "solver/TimeControl.hpp"

#include <json/json.h>
#include <range/v3/view.hpp>

#include "CubeJson.hpp"

using namespace neon;

std::string const time_json("{\"Period\" : 1.0, \"Increments\": { "
                            "\"Initial\" : 1.0, \"Minimum\" : 0.001, \"Maximum\" : 10.0 }}");

std::string const nonlinear_options_json("{\"NonlinearOptions\" : { "
                                         "\"DisplacementIncrementTolerance\" : 1.0e-5, "
                                         "\"ResidualForceTolerance\" : 0.001}}");

std::string const nonlinear_options_disp_broken_json("{\"NonlinearOptions\" : { "
                                                     "\"DisacementIncrementTolerance\" : 1.0e-5, "
                                                     "\"ResidualForceTolerance\" : 0.001}}");

std::string const nonlinear_options_force_broken_json("{\"NonlinearOptions\" : { "
                                                      "\"DisplacementIncrementTolerance\" : "
                                                      "1.0e-5, "
                                                      "\"ResidForceTolerance\" : 0.001}}");

std::string const visualisation_json("{\"Fields\" : [\"Displacement\", \"CauchyStress\"]}");

TEST_CASE("Nonlinear system equilibrium solver test")
{
    using solid::MaterialCoordinates;
    using solid::femMesh;
    using solid::femStaticMatrix;

    // Read in a cube mesh from the json input file and use this to
    // test the functionality of the basic mesh
    Json::Value mesh_data, material_data, simulation_data, solver_data, time_data,
        visualisation_data, nonlinear_data;
    Json::Reader reader;

    REQUIRE(reader.parse(cube_mesh_json().c_str(), mesh_data));
    REQUIRE(reader.parse(material_data_json().c_str(), material_data));
    REQUIRE(reader.parse(simulation_data_json().c_str(), simulation_data));
    REQUIRE(reader.parse(solver_data_json().c_str(), solver_data));
    REQUIRE(reader.parse(time_json.c_str(), time_data));
    REQUIRE(reader.parse(visualisation_json.c_str(), visualisation_data));

    // Create the test objects
    BasicMesh basic_mesh(mesh_data);
    NodalCoordinates nodal_coordinates(mesh_data);
    femMesh fem_mesh(basic_mesh, material_data, simulation_data);

    SECTION("Correct behaviour")
    {
        REQUIRE(reader.parse(nonlinear_options_json.c_str(), nonlinear_data));

        // Create the system
        femStaticMatrix fem_matrix(fem_mesh,
                                   Visualisation("UnitTest", fem_mesh, visualisation_data),
                                   solver_data,
                                   nonlinear_data["NonlinearOptions"],
                                   time_data);

        fem_matrix.solve();
    }
    SECTION("NonlinearOptions displacement increment incorrect")
    {
        REQUIRE(reader.parse(nonlinear_options_disp_broken_json.c_str(), nonlinear_data));

        REQUIRE_THROWS_AS(femStaticMatrix(fem_mesh,
                                          Visualisation("UnitTest", fem_mesh, visualisation_data),
                                          solver_data,
                                          nonlinear_data["NonlinearOptions"],
                                          time_data),
                          std::runtime_error);
    }
    SECTION("NonlinearOptions residual force incorrect")
    {
        REQUIRE(reader.parse(nonlinear_options_force_broken_json.c_str(), nonlinear_data));

        REQUIRE_THROWS_AS(femStaticMatrix(fem_mesh,
                                          Visualisation("UnitTest", fem_mesh, visualisation_data),
                                          solver_data,
                                          nonlinear_data["NonlinearOptions"],
                                          time_data),
                          std::runtime_error);
    }
}
