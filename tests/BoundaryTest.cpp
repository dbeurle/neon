
#define CATCH_CONFIG_MAIN

#include <catch.hpp>

#include "mesh/DofAllocator.hpp"

#include "interpolations/Triangle3.hpp"
#include "quadrature/TriangleQuadrature.hpp"

#include "mesh/solid/boundary/Boundary.hpp"
#include "mesh/solid/boundary/NonFollowerLoad.hpp"

#include <json/json.h>
#include <range/v3/view.hpp>

#include <memory>

using namespace neon;
using namespace ranges;

TEST_CASE("Dof List Allocation", "[DofAllocator]")
{
    SECTION("One element")
    {
        std::vector<List> const nodal_connectivity = {{0, 1, 2, 3}};

        std::vector<List> const known_dof_list{{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11}};

        std::vector<List> const computed_list = allocate_dof_list(3, nodal_connectivity);

        REQUIRE(view::set_difference(computed_list.at(0), known_dof_list.at(0)).empty());
    }
    SECTION("Two elements")
    {
        std::vector<List> const nodal_connectivity = {{0, 1, 2, 3}, {4, 2, 1, 5}};

        std::vector<List> const known_dof_list{{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11},
                                               {12, 13, 14, 6, 7, 8, 3, 4, 5, 15, 16, 17}};

        std::vector<List> const computed_list = allocate_dof_list(3, nodal_connectivity);

        REQUIRE(view::set_difference(computed_list.at(0), known_dof_list.at(0)).empty());
        REQUIRE(view::set_difference(computed_list.at(1), known_dof_list.at(1)).empty());
    }
}
TEST_CASE("Dof List Filter", "[DofAllocator]")
{
    SECTION("One element 0 offset")
    {
        std::vector<List> const nodal_connectivity = {{0, 1, 2, 3}};

        std::vector<List> const known_dof_list{{0, 3, 6, 9}};

        std::vector<List> const computed_list = filter_dof_list(3, 0, nodal_connectivity);

        REQUIRE(view::set_difference(computed_list.at(0), known_dof_list.at(0)).empty());
    }
    SECTION("One element 1 offset")
    {
        std::vector<List> const nodal_connectivity = {{0, 1, 2, 3}};

        std::vector<List> const known_dof_list{{1, 4, 7, 10}};

        std::vector<List> const computed_list = filter_dof_list(3, 1, nodal_connectivity);

        REQUIRE(view::set_difference(computed_list.at(0), known_dof_list.at(0)).empty());
    }
}
TEST_CASE("Boundary unit test", "[Boundary]")
{
    SECTION("Ramped load")
    {
        Boundary boundary(2.0, true);
        REQUIRE(boundary.interpolate_prescribed_value(0.5) == Approx(1.0));
        REQUIRE(boundary.interpolate_prescribed_value(1.0) == Approx(2.0));
    }
    SECTION("Instantaneous load")
    {
        Boundary boundary(2.0, false);
        REQUIRE(boundary.interpolate_prescribed_value(0.5) == Approx(2.0));
        REQUIRE(boundary.interpolate_prescribed_value(1.0) == Approx(2.0));
    }
}
TEST_CASE("Traction test for triangle", "[Traction]")
{
    using namespace neon::solid;

    // Build a right angled triangle
    Vector coordinates(9);
    coordinates << 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0;
    MaterialCoordinates material_coordinates(coordinates);

    std::vector<List> nodal_connectivity = {{0, 1, 2}};
    std::vector<List> dof_list = {{0, 3, 6}};

    Traction traction(nodal_connectivity,
                      std::make_shared<MaterialCoordinates>(material_coordinates),
                      1.0,  // Prescribed load
                      true, // Is load ramped
                      0,    // Dof offset
                      std::make_unique<Triangle3>(TriangleQuadrature::Rule::OnePoint));

    REQUIRE(traction.elements() == 1);

    auto const & [ dofs, t ] = traction.external_force(0, 1.0);

    REQUIRE((t - 1.0 / 6.0 * Vector3::Ones()).norm() == Approx(0.0));
    REQUIRE(view::set_difference(dof_list.at(0), dofs).empty());
}
