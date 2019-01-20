
#include <catch2/catch.hpp>

#include "interpolations/pyramid.hpp"
#include "quadrature/pyramid/bedrosian.hpp"

#include "interpolations/interpolation_factory.hpp"
#include "mesh/element_topology.hpp"

#include "io/json.hpp"

#include <numeric>

using namespace neon;

json full() { return json::parse("{\"element_options\" : {\"quadrature\" : \"full\"}}"); }
json reduced() { return json::parse("{\"element_options\" : {\"quadrature\" : \"reduced\"}}"); }

constexpr auto ZERO_MARGIN = 1.0e-5;

TEST_CASE("pyramid")
{
    using namespace neon::quadrature::pyramid;

    SECTION("Bedrosian scheme")
    {
        // Check 1 and 6 point rule
        bedrosian deg1(1);
        bedrosian deg2(2);
        bedrosian deg5(5);

        REQUIRE(deg1.points() == 1);
        REQUIRE(deg2.points() == 8);
        REQUIRE(deg5.points() == 27);

        REQUIRE(deg1.degree() >= 1);
        REQUIRE(deg2.degree() >= 2);
        REQUIRE(deg5.degree() >= 3);

        REQUIRE(std::accumulate(begin(deg1.weights()), end(deg1.weights()), 0.0) == Approx(4.0 / 3.0));
        REQUIRE(std::accumulate(begin(deg2.weights()), end(deg2.weights()), 0.0) == Approx(4.0 / 3.0));
        REQUIRE(std::accumulate(begin(deg5.weights()), end(deg5.weights()), 0.0) == Approx(4.0 / 3.0));
    }
    SECTION("Five node evaluation")
    {
        pyramid5 element;
        bedrosian scheme(1);

        REQUIRE(element.number_of_nodes() == 5);

        for (auto const& coordinate : scheme.coordinates())
        {
            auto const [N, dN] = element.evaluate(coordinate);

            REQUIRE(N.sum() == Approx(1.0));

            REQUIRE(dN.col(0).sum() == Approx(0.0).margin(ZERO_MARGIN));
            REQUIRE(dN.col(1).sum() == Approx(0.0).margin(ZERO_MARGIN));
            REQUIRE(dN.col(2).sum() == Approx(0.0).margin(ZERO_MARGIN));
        }
    }
    SECTION("Thirteen node evaluation")
    {
        pyramid13 element;
        bedrosian scheme(3);

        REQUIRE(element.number_of_nodes() == 13);

        for (auto const& coordinate : scheme.coordinates())
        {
            auto const [N, dN] = element.evaluate(coordinate);

            REQUIRE(N.sum() == Approx(1.0));

            REQUIRE(dN.col(0).sum() == Approx(0.0).margin(ZERO_MARGIN));
            REQUIRE(dN.col(1).sum() == Approx(0.0).margin(ZERO_MARGIN));
            REQUIRE(dN.col(2).sum() == Approx(0.0).margin(ZERO_MARGIN));
        }
        // REQUIRE(element.local_quadrature_extrapolation().rows() == 13);
        // REQUIRE(element.local_quadrature_extrapolation().cols() == 8);
        // REQUIRE(element.local_quadrature_extrapolation().allFinite());
    }
    // SECTION("Five node volume - one point")
    // {
    //     pyramid5 element(pyramid_quadrature::point::one);
    //
    //     auto const base = 1.0;
    //     auto const height = 1.0;
    //
    //     matrix3x x(3, 5);
    //     x << 0.0, base, base, 0.0, height / 2.0, //
    //         0.0, 0.0, base, base, height / 2.0,  //
    //         0.0, 0.0, 0.0, 0.0, height;
    //
    //     REQUIRE(element.compute_measure(x) == Approx(1.0 / 3.0));
    // }
    // SECTION("Five node volume - eight point")
    // {
    //     pyramid5 element(pyramid_quadrature::point::eight);
    //
    //     auto constexpr base = 1.0;
    //     auto constexpr height = 1.0;
    //
    //     matrix3x x(3, 5);
    //     x << 0.0, base, base, 0.0, height / 2.0, //
    //         0.0, 0.0, base, base, height / 2.0,  //
    //         0.0, 0.0, 0.0, 0.0, height;
    //
    //     REQUIRE(element.compute_measure(x) == Approx(1.0 / 3.0));
    // }
    // SECTION("Thirteen node volume - eight point")
    // {
    //     pyramid13 element(pyramid_quadrature::point::eight);
    //
    //     matrix3x x(3, 13);
    //     x << 1.0, 0.0, -1.0, -1.0, -1.0, 0.0, 1.0, 1.0, 0.5, -0.5, -0.5, 0.5, 0.0, //
    //         1.0, 1.0, 1.0, 0.0, -1.0, -1.0, -1.0, 0.0, 0.5, 0.5, -0.5, -0.5, 0.0,  //
    //         0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.5, 0.5, 0.5, 0.5, 1.0;
    //
    //     REQUIRE(element.compute_measure(x) == Approx(4.0 / 3.0));
    // }
}
