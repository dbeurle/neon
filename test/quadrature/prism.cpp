
#include <catch2/catch.hpp>

#include "quadrature/prism/felippa.hpp"
#include "quadrature/prism/kubatko.hpp"

#include "interpolations/prism.hpp"

#include "interpolations/interpolation_factory.hpp"
#include "mesh/element_topology.hpp"

#include "io/json.hpp"

#include <numeric>

using namespace neon;

json full() { return json::parse("{\"element_options\" : {\"quadrature\" : \"full\"}}"); }
json reduced() { return json::parse("{\"element_options\" : {\"quadrature\" : \"reduced\"}}"); }

constexpr auto ZERO_MARGIN = 1.0e-5;

matrix3x six_node_coordinates()
{
    matrix3x x(3, 6);
    x << 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, //
        0.0, 1.0, 0.0, 0.0, 1.0, 0.0,  //
        -1.0, -1.0, -1.0, 1.0, 1.0, 1.0;
    return x;
}

matrix3x fifteen_node_coordinates()
{
    matrix3x x(3, 15);
    x << 1.0, 0.0, 0.0, 0.5, 0.0, 0.5, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.5, 0.0, 0.5, //
        0.0, 1.0, 0.0, 0.5, 0.5, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.5, 0.5, 0.0,  //
        -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0;
    return x;
}

TEST_CASE("Prism quadrature scheme test")
{
    using namespace neon::quadrature::prism;

    SECTION("Felippa quadrature")
    {
        // Check one to six degree rules
        felippa p1(1);
        felippa p2(2);
        felippa p3(3);
        felippa p4(4);
        felippa p5(5);
        felippa p6(6);

        REQUIRE(p1.points() == 1);
        REQUIRE(p2.points() == 6);
        REQUIRE(p3.points() == 6);
        REQUIRE(p4.points() == 18);
        REQUIRE(p5.points() == 21);
        REQUIRE(p6.points() == 48);

        REQUIRE(p1.degree() == 1);
        REQUIRE(p2.degree() == 3);
        REQUIRE(p3.degree() == 3);
        REQUIRE(p4.degree() == 4);
        REQUIRE(p5.degree() == 5);
        REQUIRE(p6.degree() == 6);

        REQUIRE(std::accumulate(begin(p1.weights()), end(p1.weights()), 0.0) == Approx(1.0));
        REQUIRE(std::accumulate(begin(p2.weights()), end(p2.weights()), 0.0) == Approx(1.0));
        REQUIRE(std::accumulate(begin(p3.weights()), end(p3.weights()), 0.0) == Approx(1.0));
        REQUIRE(std::accumulate(begin(p4.weights()), end(p4.weights()), 0.0) == Approx(1.0));
        REQUIRE(std::accumulate(begin(p5.weights()), end(p5.weights()), 0.0) == Approx(1.0));
        REQUIRE(std::accumulate(begin(p6.weights()), end(p6.weights()), 0.0) == Approx(1.0));
    }
    SECTION("Kubatko quadrature")
    {
        // Check one to five degree rules
        kubatko p1(1);
        kubatko p2(2);
        kubatko p3(3);
        kubatko p4(4);
        kubatko p5(5);

        REQUIRE(p1.points() == 1);
        REQUIRE(p2.points() == 4);
        REQUIRE(p3.points() == 6);
        REQUIRE(p4.points() == 11);
        REQUIRE(p5.points() == 15);

        REQUIRE(p1.degree() >= 1);
        REQUIRE(p2.degree() >= 2);
        REQUIRE(p3.degree() >= 3);
        REQUIRE(p4.degree() >= 4);
        REQUIRE(p5.degree() >= 5);

        REQUIRE(std::accumulate(begin(p1.weights()), end(p1.weights()), 0.0) == Approx(1.0));
        REQUIRE(std::accumulate(begin(p2.weights()), end(p2.weights()), 0.0) == Approx(1.0));
        REQUIRE(std::accumulate(begin(p3.weights()), end(p3.weights()), 0.0) == Approx(1.0));
        REQUIRE(std::accumulate(begin(p4.weights()), end(p4.weights()), 0.0) == Approx(1.0));
        REQUIRE(std::accumulate(begin(p5.weights()), end(p5.weights()), 0.0) == Approx(1.0));
    }
}
TEST_CASE("Element checks")
{
    using namespace neon::quadrature::prism;

    SECTION("prism6 kubatko")
    {
        prism6 element;
        kubatko scheme(1);

        REQUIRE(element.number_of_nodes() == 6);

        for (auto const& coordinate : scheme.coordinates())
        {
            auto const [N, dN] = element.evaluate(coordinate);

            REQUIRE(N.sum() == Approx(1.0));

            REQUIRE(dN.col(0).sum() == Approx(0.0).margin(ZERO_MARGIN));
            REQUIRE(dN.col(1).sum() == Approx(0.0).margin(ZERO_MARGIN));
            REQUIRE(dN.col(2).sum() == Approx(0.0).margin(ZERO_MARGIN));
        }
    }
    SECTION("prism6 felippa")
    {
        prism6 element;
        felippa scheme(1);

        REQUIRE(element.number_of_nodes() == 6);

        for (auto const& coordinate : scheme.coordinates())
        {
            auto const [N, dN] = element.evaluate(coordinate);

            REQUIRE(N.sum() == Approx(1.0));

            REQUIRE(dN.col(0).sum() == Approx(0.0).margin(ZERO_MARGIN));
            REQUIRE(dN.col(1).sum() == Approx(0.0).margin(ZERO_MARGIN));
            REQUIRE(dN.col(2).sum() == Approx(0.0).margin(ZERO_MARGIN));
        }
    }
    // SECTION("Six node - six point area")
    // {
    //     prism6 pri6(prism_quadrature::point::six);
    //
    //     REQUIRE(pri6.compute_measure(six_node_coordinates()) == Approx(1.0));
    // }
    // SECTION("Fifteen node - six point evaluation")
    // {
    //     prism15 element(prism_quadrature::point::six);
    //
    //     REQUIRE(element.number_of_nodes() == 15);
    //     REQUIRE(element.quadrature().points() == 6);
    //
    //     element.quadrature().for_each([&](auto const& femval, auto) {
    //         auto const [N, dN] = femval;
    //
    //         REQUIRE(N.sum() == Approx(1.0));
    //
    //         REQUIRE(dN.col(0).sum() == Approx(0.0).margin(ZERO_MARGIN));
    //         REQUIRE(dN.col(1).sum() == Approx(0.0).margin(ZERO_MARGIN));
    //         REQUIRE(dN.col(2).sum() == Approx(0.0).margin(ZERO_MARGIN));
    //     });
    //     REQUIRE(element.local_quadrature_extrapolation().rows() == 15);
    //     REQUIRE(element.local_quadrature_extrapolation().cols() == 6);
    //     REQUIRE(element.local_quadrature_extrapolation().allFinite());
    // }
    // SECTION("Fifteen node - nine point area")
    // {
    //     prism15 pri15(prism_quadrature::point::nine);
    //
    //     REQUIRE(pri15.compute_measure(fifteen_node_coordinates()) == Approx(1.0));
    // }
}
