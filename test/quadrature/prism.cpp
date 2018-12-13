
#include <catch2/catch.hpp>

#include "quadrature/prism/prism_quadrature.hpp"
#include "quadrature/prism/felippa_prism.hpp"
#include "quadrature/prism/kubatko_prism.hpp"

#include "interpolations/prism.hpp"

#include "interpolations/interpolation_factory.hpp"
#include "mesh/element_topology.hpp"

#include "io/json.hpp"

#include <range/v3/numeric/accumulate.hpp>

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

TEST_CASE("Prism quadrature scheme test", "[prism_quadrature]")
{
    SECTION("Prism Gauss Quadrature")
    {
        // Check 1 and 6 point rule
        prism_quadrature p1(prism_quadrature::point::one);
        prism_quadrature p6(prism_quadrature::point::six);
        prism_quadrature p9(prism_quadrature::point::nine);

        REQUIRE(p1.points() == 1);
        REQUIRE(p6.points() == 6);
        REQUIRE(p9.points() == 9);

        REQUIRE(ranges::accumulate(p1.weights(), 0.0) == Approx(1.0));
        REQUIRE(ranges::accumulate(p6.weights(), 0.0) == Approx(1.0));
        REQUIRE(ranges::accumulate(p9.weights(), 0.0) == Approx(1.0));
    }
    SECTION("Felippa quadrature")
    {
        // Check one to six degree rules
        felippa_prism p1(1);
        felippa_prism p2(2);
        felippa_prism p3(3);
        felippa_prism p4(4);
        felippa_prism p5(5);
        felippa_prism p6(6);

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

        REQUIRE(ranges::accumulate(p1.weights(), 0.0) == Approx(1.0));
        REQUIRE(ranges::accumulate(p2.weights(), 0.0) == Approx(1.0));
        REQUIRE(ranges::accumulate(p3.weights(), 0.0) == Approx(1.0));
        REQUIRE(ranges::accumulate(p4.weights(), 0.0) == Approx(1.0));
        REQUIRE(ranges::accumulate(p5.weights(), 0.0) == Approx(1.0));
        REQUIRE(ranges::accumulate(p6.weights(), 0.0) == Approx(1.0));
    }
    SECTION("Kubatko quadrature")
    {
        // Check one to five degree rules
        kubatko_prism p1(1);
        kubatko_prism p2(2);
        kubatko_prism p3(3);
        kubatko_prism p4(4);
        kubatko_prism p5(5);

        REQUIRE(p1.points() == 1);
        REQUIRE(p2.points() == 4);
        REQUIRE(p3.points() == 6);
        REQUIRE(p4.points() == 11);
        REQUIRE(p5.points() == 15);

        REQUIRE(p1.degree() == 1);
        REQUIRE(p2.degree() == 2);
        REQUIRE(p3.degree() == 3);
        REQUIRE(p4.degree() == 4);
        REQUIRE(p5.degree() == 5);

        REQUIRE(ranges::accumulate(p1.weights(), 0.0) == Approx(1.0));
        REQUIRE(ranges::accumulate(p2.weights(), 0.0) == Approx(1.0));
        REQUIRE(ranges::accumulate(p3.weights(), 0.0) == Approx(1.0));
        REQUIRE(ranges::accumulate(p4.weights(), 0.0) == Approx(1.0));
        REQUIRE(ranges::accumulate(p5.weights(), 0.0) == Approx(1.0));
    }
}
TEST_CASE("Element checks")
{
    SECTION("Six node - one point evaluation")
    {
        prism6 element(prism_quadrature::point::one);

        REQUIRE(element.number_of_nodes() == 6);
        REQUIRE(element.quadrature().points() == 1);

        element.quadrature().for_each([&](auto const& femval, auto) {
            auto const& [N, rhea] = femval;

            REQUIRE(N.sum() == Approx(1.0));

            REQUIRE(rhea.col(0).sum() == Approx(0.0).margin(ZERO_MARGIN));
            REQUIRE(rhea.col(1).sum() == Approx(0.0).margin(ZERO_MARGIN));
            REQUIRE(rhea.col(2).sum() == Approx(0.0).margin(ZERO_MARGIN));
        });
        REQUIRE(element.local_quadrature_extrapolation().rows() == 6);
        REQUIRE(element.local_quadrature_extrapolation().cols() == 1);
        REQUIRE(element.local_quadrature_extrapolation().allFinite());
    }
    SECTION("Six node - six point evaluation")
    {
        prism6 element(prism_quadrature::point::six);

        REQUIRE(element.number_of_nodes() == 6);
        REQUIRE(element.quadrature().points() == 6);

        element.quadrature().for_each([&](auto const& femval, auto) {
            auto const& [N, rhea] = femval;

            REQUIRE(N.sum() == Approx(1.0));

            REQUIRE(rhea.col(0).sum() == Approx(0.0).margin(ZERO_MARGIN));
            REQUIRE(rhea.col(1).sum() == Approx(0.0).margin(ZERO_MARGIN));
            REQUIRE(rhea.col(2).sum() == Approx(0.0).margin(ZERO_MARGIN));
        });
        REQUIRE(element.local_quadrature_extrapolation().rows() == 6);
        REQUIRE(element.local_quadrature_extrapolation().cols() == 6);
        REQUIRE(element.local_quadrature_extrapolation().allFinite());
    }
    SECTION("Six node - six point area")
    {
        prism6 pri6(prism_quadrature::point::six);

        REQUIRE(pri6.compute_measure(six_node_coordinates()) == Approx(1.0));
    }
    SECTION("Fifteen node - six point evaluation")
    {
        prism15 element(prism_quadrature::point::six);

        REQUIRE(element.number_of_nodes() == 15);
        REQUIRE(element.quadrature().points() == 6);

        element.quadrature().for_each([&](auto const& femval, auto) {
            auto const& [N, rhea] = femval;

            REQUIRE(N.sum() == Approx(1.0));

            REQUIRE(rhea.col(0).sum() == Approx(0.0).margin(ZERO_MARGIN));
            REQUIRE(rhea.col(1).sum() == Approx(0.0).margin(ZERO_MARGIN));
            REQUIRE(rhea.col(2).sum() == Approx(0.0).margin(ZERO_MARGIN));
        });
        REQUIRE(element.local_quadrature_extrapolation().rows() == 15);
        REQUIRE(element.local_quadrature_extrapolation().cols() == 6);
        REQUIRE(element.local_quadrature_extrapolation().allFinite());
    }
    SECTION("Fifteen node - nine point area")
    {
        prism15 pri15(prism_quadrature::point::nine);

        REQUIRE(pri15.compute_measure(fifteen_node_coordinates()) == Approx(1.0));
    }
    SECTION("Virtual methods check")
    {
        REQUIRE(make_volume_interpolation(element_topology::prism6, full())->number_of_nodes() == 6);
        REQUIRE(make_volume_interpolation(element_topology::prism15, full())->number_of_nodes() == 15);
    }
}
