
#define CATCH_CONFIG_MAIN

#include <catch.hpp>

#include "mesh/DofAllocator.hpp"

#include "interpolations/Triangle3.hpp"
#include "quadrature/TriangleQuadrature.hpp"

#include "mesh/Submesh.hpp"

#include "mesh/generic/Boundary.hpp"
#include "mesh/mechanical/solid/boundary/NonFollowerLoad.hpp"

#include "mesh/diffusion/boundary/NewtonConvection.hpp"

#include <json/json.h>
#include <range/v3/view/set_algorithm.hpp>

#include <memory>

#include "CubeJson.hpp"

Json::CharReaderBuilder reader;
JSONCPP_STRING input_errors;

using namespace neon;
using namespace ranges;

constexpr auto ZERO_MARGIN = 1.0e-5;

TEST_CASE("Boundary unit test", "[Boundary]")
{
    SECTION("Check time data saved correctly")
    {
        Json::Value times, loads;

        std::istringstream times_stream("[0.0, 1.0, 2.0, 3.0]");
        std::istringstream loads_stream("[0.0, 1.0, 2.0, 3.0]");
        REQUIRE(Json::parseFromStream(reader, times_stream, &times, &input_errors));
        REQUIRE(Json::parseFromStream(reader, loads_stream, &loads, &input_errors));

        Boundary boundary(times, loads);

        auto const time_history = boundary.time_history();

        REQUIRE(time_history[0] == Approx(0.0).margin(ZERO_MARGIN));
        REQUIRE(time_history[1] == Approx(1.0));
        REQUIRE(time_history[2] == Approx(2.0));
        REQUIRE(time_history[3] == Approx(3.0));
    }
    SECTION("Monotonic loading interpolation test")
    {
        Json::Value times, loads;

        std::istringstream times_stream("[0.0, 1.0, 2.0, 3.0]");
        std::istringstream loads_stream("[0.0, 0.5, 1.0, 1.5]");
        REQUIRE(Json::parseFromStream(reader, times_stream, &times, &input_errors));
        REQUIRE(Json::parseFromStream(reader, loads_stream, &loads, &input_errors));

        Boundary boundary(times, loads);

        REQUIRE(boundary.interpolate_prescribed_load(0.75) == Approx(0.375));
        REQUIRE(boundary.interpolate_prescribed_load(0.5) == Approx(0.25));
        REQUIRE(boundary.interpolate_prescribed_load(1.0) == Approx(0.5));
        REQUIRE(boundary.interpolate_prescribed_load(1.9) == Approx(0.95));
        REQUIRE(boundary.interpolate_prescribed_load(2.0) == Approx(1.0));
        REQUIRE(boundary.interpolate_prescribed_load(2.5) == Approx(1.25));
        REQUIRE(boundary.interpolate_prescribed_load(3.0) == Approx(1.5));
        REQUIRE(boundary.interpolate_prescribed_load(2.9999999999999) == Approx(1.5));
    }
    SECTION("Unload interpolation test")
    {
        Json::Value times, loads;

        std::istringstream times_stream("[0.0, 1.0, 2.0, 3.0]");
        std::istringstream loads_stream("[0.0, 1.0, 0.0, 3.0]");
        REQUIRE(Json::parseFromStream(reader, times_stream, &times, &input_errors));
        REQUIRE(Json::parseFromStream(reader, loads_stream, &loads, &input_errors));

        Boundary boundary(times, loads);

        REQUIRE(boundary.interpolate_prescribed_load(0.0) == Approx(0.0).margin(ZERO_MARGIN));
        REQUIRE(boundary.interpolate_prescribed_load(0.5) == Approx(0.5));
        REQUIRE(boundary.interpolate_prescribed_load(1.0) == Approx(1.0));
        REQUIRE(boundary.interpolate_prescribed_load(1.5) == Approx(0.5));
        REQUIRE(boundary.interpolate_prescribed_load(2.0) == Approx(0.0).margin(ZERO_MARGIN));
    }
    SECTION("Non-matching length error test")
    {
        REQUIRE_THROWS_AS(Boundary("[0.0, 1.0, 3.0]", "[0.0, 0.5, 1.0, 1.5]"), std::runtime_error);
    }
    SECTION("Unordered time error test")
    {
        REQUIRE_THROWS_AS(Boundary("[0.0, 10.0, 3.0]", "[0.0, 0.5, 1.0]"), std::runtime_error);
    }
}
TEST_CASE("Traction test for triangle", "[Traction]")
{
    using namespace neon::mechanical::solid;

    Json::Value times, loads;

    // Build a right angled triangle
    vector coordinates(9);
    coordinates << 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0;

    auto material_coordinates = std::make_shared<MaterialCoordinates>(coordinates);

    std::vector<List> nodal_connectivity = {{0, 1, 2}};
    std::vector<List> dof_list = {{0, 3, 6}};

    SECTION("Unit load")
    {
        std::istringstream times_stream("[0.0, 1.0]");
        std::istringstream loads_stream("[0.0, 1.0]");

        REQUIRE(Json::parseFromStream(reader, times_stream, &times, &input_errors));
        REQUIRE(Json::parseFromStream(reader, loads_stream, &loads, &input_errors));

        traction patch(std::make_unique<Triangle3>(TriangleQuadrature::Rule::OnePoint),
                       nodal_connectivity,
                       dof_list,
                       material_coordinates,
                       times,
                       loads);

        REQUIRE(patch.elements() == 1);

        auto const& [dofs, t] = patch.external_force(0, 1.0);

        REQUIRE((t - 1.0 / 6.0 * vector3::Ones()).norm() == Approx(0.0).margin(ZERO_MARGIN));
        REQUIRE(view::set_difference(dof_list.at(0), dofs).empty());
    }
    SECTION("Twice unit load")
    {
        std::istringstream times_stream("[0.0, 1.0]");
        std::istringstream loads_stream("[0.0, 2.0]");
        REQUIRE(Json::parseFromStream(reader, times_stream, &times, &input_errors));
        REQUIRE(Json::parseFromStream(reader, loads_stream, &loads, &input_errors));

        traction patch(std::make_unique<Triangle3>(TriangleQuadrature::Rule::OnePoint),
                       nodal_connectivity,
                       dof_list,
                       material_coordinates,
                       times,
                       loads);

        REQUIRE(patch.elements() == 1);

        auto const& [dofs, t] = patch.external_force(0, 1.0);

        REQUIRE((t - 2.0 / 6.0 * vector3::Ones()).norm() == Approx(0.0).margin(ZERO_MARGIN));
        REQUIRE(view::set_difference(dof_list.at(0), dofs).empty());
    }
}
TEST_CASE("Pressure test for triangle", "[Pressure]")
{
    using namespace neon::mechanical::solid;

    Json::Value times, loads;

    // Build a right angled triangle
    vector coordinates(9);
    coordinates << 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0;

    auto material_coordinates = std::make_shared<MaterialCoordinates>(coordinates);

    std::vector<List> nodal_connectivity = {{0, 1, 2}};
    std::vector<List> dof_list = {{0, 1, 2, 3, 4, 5, 6, 7, 8}};

    SECTION("Unit load")
    {
        std::istringstream times_stream("[0.0, 1.0]");
        std::istringstream loads_stream("[0.0, -1.0]");

        REQUIRE(Json::parseFromStream(reader, times_stream, &times, &input_errors));
        REQUIRE(Json::parseFromStream(reader, loads_stream, &loads, &input_errors));

        pressure pressure_patch(std::make_unique<Triangle3>(TriangleQuadrature::Rule::OnePoint),
                                nodal_connectivity,
                                dof_list,
                                material_coordinates,
                                times,
                                loads);

        REQUIRE(pressure_patch.elements() == 1);

        auto const& [dofs, f_ext] = pressure_patch.external_force(0, 1.0);

        REQUIRE(view::set_difference(dof_list.at(0), dofs).empty());

        // Compare the z-component to the analytical solution
        REQUIRE((vector3(f_ext(2), f_ext(5), f_ext(8)) - 1.0 / 6.0 * vector3::Ones()).norm()
                == Approx(0.0).margin(ZERO_MARGIN));
    }
    SECTION("Twice unit load")
    {
        std::istringstream times_stream("[0.0, 1.0]");
        std::istringstream loads_stream("[0.0, -2.0]");

        REQUIRE(Json::parseFromStream(reader, times_stream, &times, &input_errors));
        REQUIRE(Json::parseFromStream(reader, loads_stream, &loads, &input_errors));

        pressure pressure_patch(std::make_unique<Triangle3>(TriangleQuadrature::Rule::OnePoint),
                                nodal_connectivity,
                                dof_list,
                                material_coordinates,
                                times,
                                loads);

        REQUIRE(pressure_patch.elements() == 1);

        auto const& [dofs, f_ext] = pressure_patch.external_force(0, 1.0);

        REQUIRE((vector3(f_ext(2), f_ext(5), f_ext(8)) - 2.0 / 6.0 * vector3::Ones()).norm()
                == Approx(0.0).margin(ZERO_MARGIN));
        REQUIRE(view::set_difference(dof_list.at(0), dofs).empty());
    }
}
TEST_CASE("Traction test for mixed mesh", "[NonFollowerLoadBoundary]")
{
    // Test the construction and population of a mixed quadrilateral and
    // triangle mesh
    using namespace neon::mechanical::solid;

    vector coordinates(3 * 5);
    coordinates << 0.0, 0.0, 0.0, //
        1.0, 0.0, 0.0,            //
        1.0, 1.0, 0.0,            //
        0.0, 1.0, 0.0,            //
        2.0, 1.0, 0.0;

    auto material_coordinates = std::make_shared<MaterialCoordinates>(coordinates);

    std::string trimesh = "{\"Name\" : \"Ysym\", \"NodalConnectivity\" : [ [ 1, 4, 2 "
                          "] ], \"Type\" : 2 }";
    std::string quadmesh = "{\"Name\" : \"Ysym\", \"NodalConnectivity\" : [ [ 0, 1, 2, 3 "
                           "] ], \"Type\" : 3 }";

    std::array<int, 3> const known_dofs_tri{{4, 13, 7}};
    std::array<int, 4> const known_dofs_quad{{1, 4, 7, 10}};

    Json::Value tri_mesh_data, quad_mesh_data, simulation_data;

    Json::Value boundary;

    std::istringstream trimesh_stream(trimesh);
    std::istringstream quadmesh_stream(quadmesh);
    std::istringstream simulation_data_stream(simulation_data_traction_json());

    std::istringstream boundary_stream("{\"Time\":[0.0, 1.0], \"Type\" : \"Traction\", "
                                       "\"Values\":{\"y\":[0.0, 1.0e-3]}}");

    REQUIRE(Json::parseFromStream(reader, trimesh_stream, &tri_mesh_data, &input_errors));
    REQUIRE(Json::parseFromStream(reader, quadmesh_stream, &quad_mesh_data, &input_errors));
    REQUIRE(Json::parseFromStream(reader, simulation_data_stream, &simulation_data, &input_errors));
    REQUIRE(Json::parseFromStream(reader, boundary_stream, &boundary, &input_errors));

    std::vector<Submesh> submeshes = {tri_mesh_data, quad_mesh_data};

    REQUIRE(submeshes.at(0).elements() == 1);
    REQUIRE(submeshes.at(1).elements() == 1);

    REQUIRE(submeshes.at(0).topology() == ElementTopology::Triangle3);
    REQUIRE(submeshes.at(1).topology() == ElementTopology::Quadrilateral4);

    REQUIRE(submeshes.at(0).nodes_per_element() == 3);
    REQUIRE(submeshes.at(1).nodes_per_element() == 4);

    // Insert this information into the nonfollower load boundary class
    // using the simulation data for the cube
    NonFollowerLoadBoundary loads(material_coordinates,
                                  submeshes,
                                  simulation_data,
                                  boundary,
                                  {{"x", 0}, {"y", 1}, {"z", 2}});

    for (auto const& [is_dof_active, meshes] : loads.interface())
    {
        if (is_dof_active)
        {
            std::visit(
                [&](auto const& mesh) {
                    auto const& [dofs_tri, f_tri] = mesh.external_force(0, 1.0);

                    REQUIRE(dofs_tri.size() == 3);
                    REQUIRE(f_tri.rows() == 3);

                    // Check dofs are correctly written
                    REQUIRE(view::set_difference(dofs_tri, known_dofs_tri).empty());
                },
                meshes.at(0));

            std::visit(
                [&](auto const& mesh) {
                    auto const& [dofs_quad, f_quad] = mesh.external_force(0, 1.0);

                    REQUIRE(dofs_quad.size() == 4);
                    REQUIRE(f_quad.rows() == 4);

                    // Check dofs are correctly written
                    REQUIRE(view::set_difference(dofs_quad, known_dofs_quad).empty());
                },
                meshes.at(1));
        }
    }
}
TEST_CASE("Newton cooling boundary conditions")
{
    using diffusion::boundary::newton_cooling;

    Json::Value times, flux, T_inf;

    // Build a right angled triangle
    vector coordinates(9);
    coordinates << 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0;

    auto material_coordinates = std::make_shared<MaterialCoordinates>(coordinates);

    std::vector<List> const nodal_connectivity = {{0, 1, 2}};
    std::vector<List> const dof_list = {{0, 3, 6}};

    SECTION("Unit load")
    {
        std::istringstream times_stream("[0.0, 1.0]"), h_stream("[0.0, 1.0]"),
            T_inf_stream("[300.0, 300.0]");

        REQUIRE(Json::parseFromStream(reader, times_stream, &times, &input_errors));
        REQUIRE(Json::parseFromStream(reader, h_stream, &flux, &input_errors));
        REQUIRE(Json::parseFromStream(reader, T_inf_stream, &T_inf, &input_errors));

        newton_cooling patch(std::make_unique<Triangle3>(TriangleQuadrature::Rule::OnePoint),
                             nodal_connectivity,
                             dof_list,
                             material_coordinates,
                             times,
                             flux,
                             T_inf);

        REQUIRE(patch.elements() == 1);

        auto const& [dofs_0, t] = patch.external_force(0, 1.0);
        auto const& [dofs_1, k] = patch.external_stiffness(0, 1.0);

        REQUIRE((t - 1.0 / 6.0 * vector3::Ones()).norm() == Approx(0.0).margin(ZERO_MARGIN));

        REQUIRE(k.norm() > 0.0);

        REQUIRE(view::set_difference(dof_list.at(0), dofs_0).empty());
        REQUIRE(view::set_difference(dof_list.at(0), dofs_1).empty());
    }
}
