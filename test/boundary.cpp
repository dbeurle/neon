
#include <catch2/catch.hpp>

#include "interpolations/triangle.hpp"
#include "quadrature/triangle/cowper.hpp"

#include "math/linear_interpolation.hpp"

#include "mesh/basic_submesh.hpp"
#include "mesh/boundary/boundary_condition.hpp"
#include "mesh/mechanics/solid/boundary/nonfollower_load.hpp"
#include "mesh/diffusion/heat/boundary/newton_convection.hpp"
#include "io/json.hpp"

#include "fixtures/cube_mesh.hpp"

#include <memory>

using namespace neon;

constexpr auto ZERO_MARGIN = 1.0e-5;

matrix3 unit_triangle_coordinates()
{
    matrix3x coordinates(3, 3);
    coordinates << 0.0, 1.0, 0.0, //
        0.0, 0.0, 1.0,            //
        0.0, 0.0, 0.0;
    return coordinates;
}

matrix3 tri_quad_coordinates()
{
    matrix3x coordinates(3, 5);
    coordinates << 0.0, 1.0, 1.0, 0.0, 2.0, //
        0.0, 0.0, 1.0, 1.0, 1.0,            //
        0.0, 0.0, 0.0, 0.0, 0.0;            //
    return coordinates;
}

TEST_CASE("linear interpolation")
{
    REQUIRE(relative_distance(2.0, 0.0, 4.0) == Approx(0.5));
    REQUIRE(relative_distance(5.0, 4.0, 6.0) == Approx(0.5));
    REQUIRE(relative_distance(-5.0, 0.0, -10.0) == Approx(0.5));

    REQUIRE(linear_interpolation(0.5, 5.0, 15.0) == Approx(10.0));
}
TEST_CASE("boundary unit test")
{
    std::vector<double> times{0.0, 1.0, 2.0, 3.0};
    SECTION("Check time data saved correctly")
    {
        std::vector<double> loads{0.0, 1.0, 2.0, 3.0};

        boundary_condition boundary(times, loads);

        auto const& time_history = boundary.times();

        REQUIRE(time_history[0] == Approx(0.0).margin(ZERO_MARGIN));
        REQUIRE(time_history[1] == Approx(1.0));
        REQUIRE(time_history[2] == Approx(2.0));
        REQUIRE(time_history[3] == Approx(3.0));
    }
    SECTION("Monotonic loading interpolation test")
    {
        std::vector<double> loads{0.0, 0.5, 1.0, 1.5};

        boundary_condition boundary(times, loads);

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

        boundary_condition boundary(times, loads);

        REQUIRE(boundary.interpolate_prescribed_load(0.0) == Approx(0.0).margin(ZERO_MARGIN));
        REQUIRE(boundary.interpolate_prescribed_load(0.5) == Approx(0.5));
        REQUIRE(boundary.interpolate_prescribed_load(1.0) == Approx(1.0));
        REQUIRE(boundary.interpolate_prescribed_load(1.5) == Approx(0.5));
        REQUIRE(boundary.interpolate_prescribed_load(2.0) == Approx(0.0).margin(ZERO_MARGIN));
    }
    SECTION("Non-matching length error test")
    {
        REQUIRE_THROWS_AS(boundary_condition(json::parse("[0.0, 1.0, 3.0]"),
                                             json::parse("[0.0, 0.5, 1.0, 1.5]")),
                          std::domain_error);
    }
    SECTION("Unordered time error test")
    {
        REQUIRE_THROWS_AS(boundary_condition(json::parse("[0.0, 10.0, 3.0]"),
                                             json::parse("[0.0, 0.5, 1.0]")),
                          std::domain_error);
    }
}
TEST_CASE("Traction test for triangle", "[Traction]")
{
    using namespace neon::mechanics::solid;

    auto coordinates = std::make_shared<material_coordinates>(unit_triangle_coordinates());

    indices nodal_connectivity(3, 1);
    nodal_connectivity << 0, 1, 2;

    indices dof_indices(3, 1);
    dof_indices << 0, 3, 6;

    SECTION("Unit load")
    {
        traction patch(nodal_connectivity,
                       dof_indices,
                       coordinates,
                       json::parse("[0.0, 1.0]"),
                       json::parse("[0.0, 1.0]"),
                       element_topology::triangle3,
                       {});

        REQUIRE(patch.elements() == 1);

        auto const dof_view = patch.local_dof_view(0);

        auto const& t = patch.external_force(0, 1.0);

        REQUIRE((t - vector3::Constant(1.0 / 6.0)).norm() == Approx(0.0).margin(ZERO_MARGIN));

        REQUIRE((dof_indices.col(0) - dof_view).sum() == 0);
    }
    SECTION("Twice unit load")
    {
        traction patch(nodal_connectivity,
                       dof_indices,
                       coordinates,
                       json::parse("[0.0, 1.0]"),
                       json::parse("[0.0, 2.0]"),
                       element_topology::triangle3,
                       {});

        REQUIRE(patch.elements() == 1);

        auto const& t = patch.external_force(0, 1.0);

        REQUIRE((t - vector3::Constant(2.0 / 6.0)).norm() == Approx(0.0).margin(ZERO_MARGIN));
        REQUIRE((dof_indices.col(0) - patch.local_dof_view(0)).sum() == 0);
    }
}
TEST_CASE("Pressure test for triangle")
{
    using namespace neon::mechanics::solid;

    auto coordinates = std::make_shared<material_coordinates>(unit_triangle_coordinates());

    indices nodal_connectivity(3, 1);
    nodal_connectivity << 0, 1, 2;

    indices dof_indices(3 * 3, 1);
    dof_indices << 0, 1, 2, 3, 4, 5, 6, 7, 8;

    SECTION("Unit load")
    {
        pressure pressure_patch(nodal_connectivity,
                                dof_indices,
                                coordinates,
                                json::parse("[0.0, 1.0]"),
                                json::parse("[0.0, -1.0]"),
                                element_topology::triangle3,
                                {});

        REQUIRE(pressure_patch.elements() == 1);

        auto const dof_view = pressure_patch.local_dof_view(0);

        auto const& f_ext = pressure_patch.external_force(0, 1.0);

        REQUIRE((dof_indices.col(0) - dof_view).sum() == 0);

        // Compare the z-component to the analytical solution
        REQUIRE((vector3(f_ext(2), f_ext(5), f_ext(8)) - vector3::Constant(1.0 / 6.0)).norm()
                == Approx(0.0).margin(ZERO_MARGIN));
    }
    SECTION("Twice unit load")
    {
        pressure pressure_patch(nodal_connectivity,
                                dof_indices,
                                coordinates,
                                json::parse("[0.0, 1.0]"),
                                json::parse("[0.0, -2.0]"),
                                element_topology::triangle3,
                                {});

        REQUIRE(pressure_patch.elements() == 1);

        auto const dof_view = pressure_patch.local_dof_view(0);

        auto const& f_ext = pressure_patch.external_force(0, 1.0);

        REQUIRE((vector3(f_ext(2), f_ext(5), f_ext(8)) - 2.0 / 6.0 * vector3::Ones()).norm()
                == Approx(0.0).margin(ZERO_MARGIN));

        REQUIRE((dof_indices.col(0) - dof_view).sum() == 0);
    }
}
TEST_CASE("Traction test for mixed mesh")
{
    // Test construction of a mixed quadrilateral and triangle mesh
    using namespace neon::mechanics::solid;

    auto coordinates = std::make_shared<material_coordinates>(tri_quad_coordinates());

    auto tri_mesh = json::parse("{\"Name\":\"Ysym\",\"NodalConnectivity\":[[1,4,2]],\"Type\":"
                                "2}");
    auto quad_mesh = json::parse("{\"Name\":\"Ysym\",\"NodalConnectivity\":[[0,1,2,3]],"
                                 "\"Type\":3}");

    indices known_dofs_tri(3, 1);
    known_dofs_tri << 4, 13, 7;

    indices known_dofs_quad(4, 1);
    known_dofs_quad << 1, 4, 7, 10;

    auto simulation_data = json::parse("{\"boundaries\" : [ "
                                       "{\"name\" : \"Ysym\", "
                                       "\"type\" : \"traction\", "
                                       "\"time\" : [0.0, 1.0],"
                                       "\"y\" : [0.0, 1.0e-3]} ], "
                                       "\"constitutive\" : {\"name\":\"neohooke\"}, "
                                       "\"element_options\" : {\"quadrature\" : \"full\"}, "
                                       "\"name\" : \"cube\"}");

    auto boundary_data = json::parse("{\"time\" : [0.0, 1.0],"
                                     "\"type\" : \"traction\","
                                     "\"y\" : [0.0, 1.0e-3]}");

    std::vector<basic_submesh> submeshes = {tri_mesh, quad_mesh};

    REQUIRE(submeshes.at(0).elements() == 1);
    REQUIRE(submeshes.at(1).elements() == 1);

    REQUIRE(submeshes.at(0).topology() == element_topology::triangle3);
    REQUIRE(submeshes.at(1).topology() == element_topology::quadrilateral4);

    REQUIRE(submeshes.at(0).nodes_per_element() == 3);
    REQUIRE(submeshes.at(1).nodes_per_element() == 4);

    // Insert this information into the nonfollower load boundary class
    // using the simulation data for the cube
    nonfollower_load_boundary loads(coordinates,
                                    submeshes,
                                    simulation_data,
                                    boundary_data,
                                    {{"x", 0}, {"y", 1}, {"z", 2}},
                                    1.0);
    // Check the triangle element
    std::visit(
        [&](auto const& mesh) {
            auto const dofs_tri = mesh.local_dof_view(0);
            auto const& f_tri = mesh.external_force(0, 1.0);

            REQUIRE(dofs_tri.size() == 3);
            REQUIRE(f_tri.rows() == 3);

            // Check dofs are correctly written
            REQUIRE((dofs_tri - known_dofs_tri.col(0)).sum() == 0);
        },
        loads.natural_interface().at(0));

    // Check the quadrilateral element
    std::visit(
        [&](auto const& mesh) {
            auto const dofs_quad = mesh.local_dof_view(0);

            auto const& f_quad = mesh.external_force(0, 1.0);

            REQUIRE(dofs_quad.size() == 4);
            REQUIRE(f_quad.rows() == 4);

            // Check dofs are correctly written
            REQUIRE((dofs_quad - known_dofs_quad.col(0)).sum() == 0);
        },
        loads.natural_interface().at(1));
}
TEST_CASE("Newton cooling boundary conditions")
{
    using diffusion::boundary::newton_convection;

    auto coordinates = std::make_shared<material_coordinates>(unit_triangle_coordinates());

    indices nodal_connectivity(3, 1);
    nodal_connectivity << 0, 1, 2;

    indices dof_indices(3, 1);
    dof_indices << 0, 3, 6;

    SECTION("Unit load")
    {
        newton_convection patch(nodal_connectivity,
                                dof_indices,
                                coordinates,
                                json::parse("[0.0, 1.0]"),
                                json::parse("[0.0, 1.0]"),
                                json::parse("[300.0, 300.0]"),
                                element_topology::triangle3,
                                {});

        REQUIRE(patch.elements() == 1);

        auto const dof_view = patch.local_dof_view(0);

        auto const& t = patch.external_force(0, 1.0);
        auto const& k = patch.external_stiffness(0, 1.0);

        REQUIRE((t - vector3::Constant(1.0 / 6.0)).norm() == Approx(0.0).margin(ZERO_MARGIN));
        REQUIRE(k.norm() > 0.0);
        REQUIRE((dof_indices - dof_view).sum() == 0);
    }
}
