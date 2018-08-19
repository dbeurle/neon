
#include <catch.hpp>

#include "mesh/nodal_coordinates.hpp"

constexpr auto ZERO_MARGIN = 1.0e-5;

TEST_CASE("Nodal coordinates")
{
    using namespace neon;

    // Build a right angled triangle (three nodes)
    matrix3x coordinates(3, 3);
    coordinates << 0.0, 1.0, 0.0, //
        0.0, 0.0, 1.0,            //
        0.0, 0.0, 0.0;

    // Setup the test case
    nodal_coordinates nodes(coordinates);

    std::vector<std::int32_t> const node_list{0, 1, 2};

    REQUIRE((nodes.coordinates() - coordinates).norm() == Approx(0.0).margin(ZERO_MARGIN));

    REQUIRE((nodes.coordinates(node_list) - coordinates).norm() == Approx(0.0).margin(ZERO_MARGIN));

    REQUIRE(nodes.coordinates(node_list).rows() == 3);
    REQUIRE(nodes.coordinates(node_list).cols() == 3);
}
