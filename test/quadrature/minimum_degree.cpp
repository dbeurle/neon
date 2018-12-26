
#include <catch2/catch.hpp>

#include "quadrature/minimum_degree.hpp"

TEST_CASE("Energy norm (mass matrix, elasticity and heat conduction)")
{
    // quadrilaterial element
    REQUIRE(neon::minimum_degree(2, 1, 1) == 1);
    // hexahedron element
    REQUIRE(neon::minimum_degree(2, 3, 1) == 3);
    // triangular element
    REQUIRE(neon::minimum_degree(2, 2, 1) == 2);
    // mass for quadratic tetrahedron
    REQUIRE(neon::minimum_degree(2, 2, 0) == 4);
    // mass for trilinear hexahedron
    REQUIRE(neon::minimum_degree(2, 3, 0) == 5);
}
