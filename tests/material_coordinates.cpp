
#define CATCH_CONFIG_MAIN

#include <catch.hpp>

#include <numeric>

#include "mesh/material_coordinates.hpp"

constexpr auto ZERO_MARGIN = 1.0e-5;

TEST_CASE("Material coordinates inheritence")
{
    using namespace neon;

    // Build a right angled triangle (three nodes)
    matrix3x coordinates(3, 3);
    coordinates << 0.0, 1.0, 0.0, //
        0.0, 0.0, 1.0,            //
        0.0, 0.0, 0.0;

    // Setup the test case
    material_coordinates single_triangle(coordinates);

    local_indices const node_list{0, 1, 2};

    REQUIRE((single_triangle.coordinates() - coordinates).norm() == Approx(0.0).margin(ZERO_MARGIN));

    REQUIRE((single_triangle.coordinates(node_list) - coordinates).norm()
            == Approx(0.0).margin(ZERO_MARGIN));

    REQUIRE(single_triangle.coordinates(node_list).rows() == 3);
    REQUIRE(single_triangle.coordinates(node_list).cols() == 3);
}
TEST_CASE("Updating material coordinates")
{
    using namespace neon;

    // Build a right angled triangle (three nodes)
    matrix3x coordinates(3, 3);
    coordinates << 0.0, 1.0, 0.0, //
        0.0, 0.0, 1.0,            //
        0.0, 0.0, 0.0;

    material_coordinates single_triangle(coordinates);

    vector test_displacements = vector::Constant(9, 1.0);

    local_indices const node_list{0, 1, 2};

    local_indices local_dof_list(9);
    std::iota(std::begin(local_dof_list), std::end(local_dof_list), 0);

    SECTION("Check initial displacements are zero")
    {
        REQUIRE(single_triangle.displacement().norm() == Approx(0.0).margin(ZERO_MARGIN));
    }
    SECTION("Test update of coordinates")
    {
        single_triangle.update_current_configuration(test_displacements);

        REQUIRE((single_triangle.displacement() - matrix3::Ones()).norm()
                == Approx(0.0).margin(ZERO_MARGIN));
    }
    SECTION("Test element initial and current configuration")
    {
        single_triangle.update_current_configuration(test_displacements);

        matrix3 known_x;
        known_x << 1, 2, 1, 1, 1, 2, 1, 1, 1;

        REQUIRE((single_triangle.initial_configuration(node_list) - coordinates).norm()
                == Approx(0.0).margin(ZERO_MARGIN));

        REQUIRE((single_triangle.current_configuration(node_list) - known_x).norm()
                == Approx(0.0).margin(ZERO_MARGIN));
    }
}
