
#include <catch.hpp>

#include "mesh/nodal_variables.hpp"

constexpr auto ZERO_MARGIN = 1.0e-5;

TEST_CASE("Scalar variable")
{
    neon::nodal_variables<1> one(10);

    REQUIRE(one.variables_per_node == 1);
    REQUIRE(one.size() == 10);

    one.update(neon::vector::Ones(10));
    REQUIRE(one.view(Eigen::placeholders::all).sum() == 10);
}

TEST_CASE("Vector 2 variable")
{
    neon::nodal_variables<2> two(10);

    REQUIRE(two.variables_per_node == 2);
    REQUIRE(two.size() == 20);

    two.update(neon::vector::Ones(20));
    REQUIRE(two.view(Eigen::placeholders::all).sum() == 20);
}

TEST_CASE("Vector 3 variable")
{
    neon::nodal_variables<3> three(10);

    REQUIRE(three.variables_per_node == 3);
    REQUIRE(three.size() == 30);

    three.update(neon::vector::Ones(30));
    REQUIRE(three.view(Eigen::placeholders::all).sum() == 30);
}

TEST_CASE("Vector 6 variable")
{
    neon::nodal_variables<6> six(10);

    REQUIRE(six.variables_per_node == 6);
    REQUIRE(six.size() == 60);

    six.update(neon::vector::Ones(60));
    REQUIRE(six.view(Eigen::placeholders::all).sum() == 60);
}
