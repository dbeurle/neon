
#include <catch2/catch.hpp>

#include "quadrature/triangle/cowper.hpp"
#include "interpolations/triangle.hpp"

#include "mesh/element_topology.hpp"
#include "interpolations/interpolation_factory.hpp"

#include "io/json.hpp"

#include <numeric>

constexpr auto ZERO_MARGIN = 1.0e-5;

using namespace neon;

json full() { return json::parse("{\"element_options\" : {\"quadrature\" : \"full\"}}"); }
json reduced() { return json::parse("{\"element_options\" : {\"quadrature\" : \"reduced\"}}"); }

TEST_CASE("Triangle quadrature")
{
    using namespace neon::quadrature::triangle;

    SECTION("triangle")
    {
        cowper deg1(1);
        cowper deg3(3);
        cowper deg5(5);

        REQUIRE(deg1.points() == 1);
        REQUIRE(deg3.points() == 3);
        REQUIRE(deg5.points() == 4);

        REQUIRE(deg1.degree() >= 1);
        REQUIRE(deg3.degree() >= 3);
        REQUIRE(deg5.degree() >= 4);

        REQUIRE(std::accumulate(begin(deg1.weights()), end(deg1.weights()), 0.0) == Approx(0.5));
        REQUIRE(std::accumulate(begin(deg3.weights()), end(deg3.weights()), 0.0) == Approx(0.5));
        REQUIRE(std::accumulate(begin(deg5.weights()), end(deg5.weights()), 0.0) == Approx(0.5));
    }
    SECTION("triangle3 cowper")
    {
        triangle3 tri3;

        cowper scheme(1);

        REQUIRE(tri3.number_of_nodes() == 3);

        for (auto const& coordinate : scheme.coordinates())
        {
            auto const [N, dN] = tri3.evaluate(coordinate);

            REQUIRE(N.sum() == Approx(1.0));

            REQUIRE(dN.col(0).sum() == Approx(0.0).margin(ZERO_MARGIN));
            REQUIRE(dN.col(1).sum() == Approx(0.0).margin(ZERO_MARGIN));
        }
    }
    SECTION("triangle6 cowper")
    {
        triangle6 tri6;
        cowper scheme(2);

        REQUIRE(tri6.number_of_nodes() == 6);

        for (auto const& coordinate : scheme.coordinates())
        {
            auto const [N, dN] = tri6.evaluate(coordinate);

            REQUIRE(N.sum() == Approx(1.0));

            REQUIRE(dN.col(0).sum() == Approx(0.0).margin(ZERO_MARGIN));
            REQUIRE(dN.col(1).sum() == Approx(0.0).margin(ZERO_MARGIN));
        }
    }
}

// SECTION("triangle3 surface area")
// {
//     triangle3 tri3(triangle_quadrature::point::one);
//
//     matrix x(3, 3);
//     x << 0.0, 0.0, 0.0, //
//         1.0, 0.0, 0.0,  //
//         0.0, 1.0, 0.0;
//     x.transposeInPlace();
//
//     REQUIRE(tri3.compute_measure(x) == Approx(0.5));
// }
// SECTION("triangle6 surface area")
// {
//     triangle6 patch(triangle_quadrature::point::three);
//
//     matrix x(6, 3);
//     x << 0.0, 0.0, 0.0, //
//         1.0, 0.0, 0.0,  //
//         0.0, 1.0, 0.0,  //
//         0.5, 0.0, 0.0,  //
//         0.5, 0.5, 0.0,  //
//         0.0, 0.5, 0.0;
//     x.transposeInPlace();
//
//     REQUIRE(patch.compute_measure(x) == Approx(0.5));
// }
