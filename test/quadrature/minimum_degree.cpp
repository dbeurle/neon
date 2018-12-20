
#include <catch2/catch.hpp>

#include "quadrature/minimum_degree.hpp"

TEST_CASE("Energy norm (elasticity or heat conduction)")
{
    // quadrilaterial element
    REQUIRE(neon::minimum_degree(2, 1, 1) == 1);
    // hexahedron element
    REQUIRE(neon::minimum_degree(2, 3, 1) == 3);
    // triangular element
    REQUIRE(neon::minimum_degree(2, 2, 1) == 2);
}
