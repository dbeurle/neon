
#define CATCH_CONFIG_MAIN

#include <catch.hpp>

#include "mesh/DofAllocator.hpp"

#include "interpolations/triangle.hpp"
#include "quadrature/triangle_quadrature.hpp"

#include "mesh/Submesh.hpp"

#include "mesh/generic/Boundary.hpp"
#include "mesh/mechanical/solid/boundary/NonFollowerLoad.hpp"

#include "mesh/diffusion/boundary/NewtonConvection.hpp"

#include "io/json.hpp"

#include <range/v3/view/set_algorithm.hpp>

#include <memory>

#include "CubeJson.hpp"

using namespace neon;
using namespace ranges;

constexpr auto ZERO_MARGIN = 1.0e-5;

TEST_CASE("Boundary unit test", "[Boundary]")
{
    std::vector<double> times{0.0, 1.0, 2.0, 3.0};
    SECTION("Check time data saved correctly")
    {
        std::vector<double> loads{0.0, 1.0, 2.0, 3.0};

        Boundary boundary(times, loads);

        auto const time_history = boundary.time_history();

        REQUIRE(time_history[0] == Approx(0.0).margin(ZERO_MARGIN));
        REQUIRE(time_history[1] == Approx(1.0));
        REQUIRE(time_history[2] == Approx(2.0));
        REQUIRE(time_history[3] == Approx(3.0));
    }
    SECTION("Monotonic loading interpolation test")
    {
        std::vector<double> loads{0.0, 0.5, 1.0, 1.5};

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
        std::vector<double> loads{0.0, 1.0, 0.0, 3.0};

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

    // Build a right angled triangle
    matrix3x coordinates(3, 3);
    coordinates << 0.0, 1.0, 0.0, //
        0.0, 0.0, 1.0,            //
        0.0, 0.0, 0.0;

    auto material_coordinates = std::make_shared<MaterialCoordinates>(coordinates);

    std::vector<List> nodal_connectivity = {{0, 1, 2}};
    std::vector<List> dof_list = {{0, 3, 6}};

    SECTION("Unit load")
    {
        traction patch(std::make_unique<triangle3>(triangle_quadrature::Rule::OnePoint),
                       nodal_connectivity,
                       dof_list,
                       material_coordinates,
                       json::parse("[0.0, 1.0]"),
                       json::parse("[0.0, 1.0]"));

        REQUIRE(patch.elements() == 1);

        auto const& [dofs, t] = patch.external_force(0, 1.0);

        REQUIRE((t - vector3::Constant(1.0 / 6.0)).norm() == Approx(0.0).margin(ZERO_MARGIN));
        REQUIRE(view::set_difference(dof_list.at(0), dofs).empty());
    }
    SECTION("Twice unit load")
    {
        traction patch(std::make_unique<triangle3>(triangle_quadrature::Rule::OnePoint),
                       nodal_connectivity,
                       dof_list,
                       material_coordinates,
                       json::parse("[0.0, 1.0]"),
                       json::parse("[0.0, 2.0]"));

        REQUIRE(patch.elements() == 1);

        auto const& [dofs, t] = patch.external_force(0, 1.0);

        REQUIRE((t - vector3::Constant(2.0 / 6.0)).norm() == Approx(0.0).margin(ZERO_MARGIN));
        REQUIRE(view::set_difference(dof_list.at(0), dofs).empty());
    }
}
TEST_CASE("Pressure test for triangle", "[Pressure]")
{
    using namespace neon::mechanical::solid;

    // Build a right angled triangle
    matrix3x coordinates(3, 3);
    coordinates << 0.0, 1.0, 0.0, //
        0.0, 0.0, 1.0,            //
        0.0, 0.0, 0.0;

    auto material_coordinates = std::make_shared<MaterialCoordinates>(coordinates);

    std::vector<List> nodal_connectivity = {{0, 1, 2}};
    std::vector<List> dof_list = {{0, 1, 2, 3, 4, 5, 6, 7, 8}};

    SECTION("Unit load")
    {
        pressure pressure_patch(std::make_unique<triangle3>(triangle_quadrature::Rule::OnePoint),
                                nodal_connectivity,
                                dof_list,
                                material_coordinates,
                                json::parse("[0.0, 1.0]"),
                                json::parse("[0.0, -1.0]"));

        REQUIRE(pressure_patch.elements() == 1);

        auto const& [dofs, f_ext] = pressure_patch.external_force(0, 1.0);

        REQUIRE(view::set_difference(dof_list.at(0), dofs).empty());

        // Compare the z-component to the analytical solution
        REQUIRE((vector3(f_ext(2), f_ext(5), f_ext(8)) - 1.0 / 6.0 * vector3::Ones()).norm()
                == Approx(0.0).margin(ZERO_MARGIN));
    }
    SECTION("Twice unit load")
    {
        pressure pressure_patch(std::make_unique<triangle3>(triangle_quadrature::Rule::OnePoint),
                                nodal_connectivity,
                                dof_list,
                                material_coordinates,
                                json::parse("[0.0, 1.0]"),
                                json::parse("[0.0, -2.0]"));

        REQUIRE(pressure_patch.elements() == 1);

        auto const& [dofs, f_ext] = pressure_patch.external_force(0, 1.0);

        REQUIRE((vector3(f_ext(2), f_ext(5), f_ext(8)) - 2.0 / 6.0 * vector3::Ones()).norm()
                == Approx(0.0).margin(ZERO_MARGIN));
        REQUIRE(view::set_difference(dof_list.at(0), dofs).empty());
    }
}
TEST_CASE("Traction test for mixed mesh", "[NonFollowerLoadBoundary]")
{
    // Test construction of a mixed quadrilateral and triangle mesh
    using namespace neon::mechanical::solid;

    matrix3x coordinates(3, 5);
    coordinates << 0.0, 1.0, 1.0, 0.0, 2.0, //
        0.0, 0.0, 1.0, 1.0, 1.0,            //
        0.0, 0.0, 0.0, 0.0, 0.0;            //

    auto material_coordinates = std::make_shared<MaterialCoordinates>(coordinates);

    auto tri_mesh = json::parse("{\"Name\":\"Ysym\",\"NodalConnectivity\":[[1,4,2]],\"Type\":"
                                "2}");
    auto quad_mesh = json::parse("{\"Name\":\"Ysym\",\"NodalConnectivity\":[[0,1,2,3]],"
                                 "\"Type\":3}");

    std::array<int, 3> const known_dofs_tri{{4, 13, 7}};
    std::array<int, 4> const known_dofs_quad{{1, 4, 7, 10}};

    auto simulation_data = json::parse("{\"BoundaryConditions\" : [ "
                                       "{\"Name\" : \"Ysym\", "
                                       "\"Type\" : \"Traction\", "
                                       "\"Time\" : [0.0, 1.0],"
                                       "\"Values\" : {\"y\" : [0.0, 1.0e-3]}} ], "
                                       "\"ConstitutiveModel\" : {\"Name\":\"NeoHooke\"}, "
                                       "\"ElementOptions\" : {\"Quadrature\" : \"Full\"}, "
                                       "\"Name\" : \"cube\"}");

    auto boundary = json::parse("{\"Time\":[0.0, "
                                "1.0],\"Type\":\"Traction\",\"Values\":{\"y\":[0.0, "
                                "1.0e-3]}}");

    std::vector<Submesh> submeshes = {tri_mesh, quad_mesh};

    REQUIRE(submeshes.at(0).elements() == 1);
    REQUIRE(submeshes.at(1).elements() == 1);

    REQUIRE(submeshes.at(0).topology() == element_topology::triangle3);
    REQUIRE(submeshes.at(1).topology() == element_topology::quadrilateral4);

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

    // Build a right angled triangle
    matrix3x coordinates(3, 3);
    coordinates << 0.0, 1.0, 0.0, //
        0.0, 0.0, 1.0,            //
        0.0, 0.0, 0.0;

    auto material_coordinates = std::make_shared<MaterialCoordinates>(coordinates);

    std::vector<List> const nodal_connectivity = {{0, 1, 2}};
    std::vector<List> const dof_list = {{0, 3, 6}};

    SECTION("Unit load")
    {
        newton_cooling patch(std::make_unique<triangle3>(triangle_quadrature::Rule::OnePoint),
                             nodal_connectivity,
                             dof_list,
                             material_coordinates,
                             json::parse("[0.0, 1.0]"),
                             json::parse("[0.0, 1.0]"),
                             json::parse("[300.0, 300.0]"));

        REQUIRE(patch.elements() == 1);

        auto const& [dofs_0, t] = patch.external_force(0, 1.0);
        auto const& [dofs_1, k] = patch.external_stiffness(0, 1.0);

        REQUIRE((t - vector3::Constant(1.0 / 6.0)).norm() == Approx(0.0).margin(ZERO_MARGIN));

        REQUIRE(k.norm() > 0.0);

        REQUIRE(view::set_difference(dof_list.at(0), dofs_0).empty());
        REQUIRE(view::set_difference(dof_list.at(0), dofs_1).empty());
    }
}
