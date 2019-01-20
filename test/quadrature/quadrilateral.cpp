
#include <catch2/catch.hpp>

#include "quadrature/quadrilateral/gauss_legendre.hpp"
#include "interpolations/quadrilateral.hpp"

#include "mesh/element_topology.hpp"
#include "interpolations/interpolation_factory.hpp"

#include <numeric>

constexpr auto ZERO_MARGIN = 1.0e-5;

using namespace neon;

TEST_CASE("Quadrilateral quadrature scheme test")
{
    using namespace neon::quadrature::quadrilateral;

    SECTION("Quadrilateral Gauss Quadrature")
    {
        gauss_legendre q1(1);
        gauss_legendre q4(3);
        gauss_legendre q9(5);

        REQUIRE(q1.points() == 1);
        REQUIRE(q4.points() == 4);
        REQUIRE(q9.points() == 9);

        REQUIRE(std::accumulate(begin(q1.weights()), end(q1.weights()), 0.0) == Approx(4.0));
        REQUIRE(std::accumulate(begin(q4.weights()), end(q4.weights()), 0.0) == Approx(4.0));
        REQUIRE(std::accumulate(begin(q9.weights()), end(q9.weights()), 0.0) == Approx(4.0));
    }
    SECTION("quadrilateral4", "[quadrilateral4.gauss_legendre]")
    {
        quadrilateral4 quad4;
        gauss_legendre integration(2);

        REQUIRE(quad4.number_of_nodes() == 4);

        for (auto const& coordinate : integration.coordinates())
        {
            auto const [N, dN] = quad4.evaluate(coordinate);

            REQUIRE(N.sum() == Approx(1.0));

            REQUIRE(dN.col(0).sum() == Approx(0.0).margin(ZERO_MARGIN));
            REQUIRE(dN.col(1).sum() == Approx(0.0).margin(ZERO_MARGIN));
        }
    }
    SECTION("quadrilateral8", "[quadrilateral8.gauss_legendre]")
    {
        quadrilateral8 quad8;
        gauss_legendre integration(2);

        REQUIRE(quad8.number_of_nodes() == 8);

        for (auto const& coordinate : integration.coordinates())
        {
            auto const [N, dN] = quad8.evaluate(coordinate);

            REQUIRE(N.sum() == Approx(1.0));

            REQUIRE(dN.col(0).sum() == Approx(0.0).margin(ZERO_MARGIN));
            REQUIRE(dN.col(1).sum() == Approx(0.0).margin(ZERO_MARGIN));
        }
    }
    SECTION("quadrilateral9", "[quadrilateral9.gauss_legendre]")
    {
        quadrilateral9 quad9;
        gauss_legendre integration(3);

        REQUIRE(quad9.number_of_nodes() == 9);

        for (auto const& coordinate : integration.coordinates())
        {
            auto const [N, dN] = quad9.evaluate(coordinate);

            REQUIRE(N.sum() == Approx(1.0));

            REQUIRE(dN.col(0).sum() == Approx(0.0).margin(ZERO_MARGIN));
            REQUIRE(dN.col(1).sum() == Approx(0.0).margin(ZERO_MARGIN));
        }
    }
    // SECTION("quadrilateral4 surface area - one point")
    // {
    //     quadrilateral4 quad4;
    //
    //     matrix x(3, 4);
    //     x << 0.0, 1.0, 1.0, 0.0, //
    //         0.0, 0.0, 1.0, 1.0,  //
    //         0.0, 0.0, 0.0, 0.0;
    //
    //     REQUIRE(quad4.compute_measure(x) == Approx(1.0));
    // }
    // SECTION("quadrilateral4 surface area - four point")
    // {
    //     quadrilateral4 quad4(quadril);
    //
    //     matrix x(3, 4);
    //     x << 0.0, 1.0, 1.0, 0.0, //
    //         0.0, 0.0, 1.0, 1.0,  //
    //         0.0, 0.0, 0.0, 0.0;
    //
    //     REQUIRE(quad4.compute_measure(x) == Approx(1.0));
    // }
    // SECTION("quadrilateral8 surface area - four point")
    // {
    //     quadrilateral8 quad8(quadril);
    //
    //     matrix x(3, 8);
    //     x << 0.0, 1.0, 1.0, 0.0, 0.5, 1.0, 0.5, 0.0, //
    //         0.0, 0.0, 1.0, 1.0, 0.0, 0.5, 1.0, 0.5,  //
    //         0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0;
    //
    //     REQUIRE(quad8.compute_measure(x) == Approx(1.0));
    // }
    // SECTION("quadrilateral8 interpolation function - nine point")
    // {
    //     quadrilateral8 quad8(quadril);
    //
    //     matrix x(3, 8);
    //     x << 0.0, 1.0, 1.0, 0.0, 0.5, 1.0, 0.5, 0.0, //
    //         0.0, 0.0, 1.0, 1.0, 0.0, 0.5, 1.0, 0.5,  //
    //         0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0;
    //
    //     REQUIRE(quad8.compute_measure(x) == Approx(1.0));
    // }
    // SECTION("quadrilateral9 surface area - four point")
    // {
    //     quadrilateral9 quad9(quadril);
    //
    //     matrix x(3, 9);
    //     x << 0.0, 1.0, 1.0, 0.0, 0.5, 1.0, 0.5, 0.0, 0.5, //
    //         0.0, 0.0, 1.0, 1.0, 0.0, 0.5, 1.0, 0.5, 0.5,  //
    //         0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0;
    //
    //     REQUIRE(quad9.compute_measure(x) == Approx(1.0));
    // }
    // SECTION("quadrilateral9 interpolation function - nine point")
    // {
    //     quadrilateral9 quad9(quadril);
    //
    //     matrix x(3, 9);
    //     x << 0.0, 1.0, 1.0, 0.0, 0.5, 1.0, 0.5, 0.0, 0.5, //
    //         0.0, 0.0, 1.0, 1.0, 0.0, 0.5, 1.0, 0.5, 0.5,  //
    //         0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0;
    //
    //     REQUIRE(quad9.compute_measure(x) == Approx(1.0));
    // }
}
