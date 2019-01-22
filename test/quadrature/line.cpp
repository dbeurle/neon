
#include <catch2/catch.hpp>

#include "quadrature/line/gauss_legendre.hpp"
#include "interpolations/line.hpp"

#include "mesh/element_topology.hpp"
#include "interpolations/interpolation_factory.hpp"

#include "io/json.hpp"

#include <numeric>

constexpr auto ZERO_MARGIN = 1.0e-5;

using namespace neon;

json full() { return json::parse("{\"element_options\" : {\"quadrature\" : \"full\"}}"); }
json reduced() { return json::parse("{\"element_options\" : {\"quadrature\" : \"reduced\"}}"); }

TEST_CASE("line quadrature scheme test")
{
    using namespace neon::quadrature::line;

    SECTION("quadrature")
    {
        // Check 1 and 3 point rule
        gauss_legendre deg1(1);
        gauss_legendre deg3(3);
        gauss_legendre deg5(5);

        REQUIRE(deg1.points() == 1);
        REQUIRE(deg3.points() == 2);
        REQUIRE(deg5.points() == 3);

        REQUIRE(deg1.degree() >= 1);
        REQUIRE(deg3.degree() >= 3);
        REQUIRE(deg5.degree() >= 5);

        REQUIRE(std::accumulate(begin(deg1.weights()), end(deg1.weights()), 0.0) == Approx(2.0));
        REQUIRE(std::accumulate(begin(deg3.weights()), end(deg3.weights()), 0.0) == Approx(2.0));
        REQUIRE(std::accumulate(begin(deg5.weights()), end(deg5.weights()), 0.0) == Approx(2.0));
    }
    SECTION("line 2 gauss legendre")
    {
        line2 line;

        gauss_legendre scheme(1);

        REQUIRE(line.number_of_nodes() == 2);

        for (auto const& coordinate : scheme.coordinates())
        {
            auto const [N, dN] = line.evaluate(coordinate);

            REQUIRE(N.sum() == Approx(1.0));
            REQUIRE(dN.sum() == Approx(0.0).margin(ZERO_MARGIN));
        }
    }
    SECTION("line 3 interpolation function - two point")
    {
        line3 line;
        gauss_legendre scheme(3);

        REQUIRE(line.number_of_nodes() == 3);

        for (auto const& coordinate : scheme.coordinates())
        {
            auto const [N, dN] = line.evaluate(coordinate);

            REQUIRE(N.sum() == Approx(1.0));

            REQUIRE(dN.sum() == Approx(0.0).margin(ZERO_MARGIN));
        }
    }
}

// SECTION("line 2 length")
// {
//     line2 line(line_quadrature::point::one);
//
//     matrix x(3, 2);
//     x << 0.0, 1.0, //
//         0.0, 0.0,  //
//         0.0, 0.0;
//
//     REQUIRE(line.compute_measure(x) == Approx(1.0));
// }
// SECTION("line 3 length")
// {
//     line3 patch(line_quadrature::point::three);
//
//     matrix x(3, 3);
//     x << 0.0, 0.5, 1.0, //
//         0.0, 0.0, 0.0,  //
//         0.0, 0.0, 0.0;
//
//     REQUIRE(patch.compute_measure(x) == Approx(1.0));
// }
