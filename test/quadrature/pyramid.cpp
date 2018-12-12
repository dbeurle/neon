
#include <catch2/catch.hpp>

#include "interpolations/pyramid.hpp"

#include "interpolations/interpolation_factory.hpp"
#include "mesh/element_topology.hpp"

#include "io/json.hpp"

#include <range/v3/numeric/accumulate.hpp>

using namespace neon;

json full() { return json::parse("{\"element_options\" : {\"quadrature\" : \"full\"}}"); }
json reduced() { return json::parse("{\"element_options\" : {\"quadrature\" : \"reduced\"}}"); }

constexpr auto ZERO_MARGIN = 1.0e-5;

TEST_CASE("Pyramid quadrature scheme test")
{
    SECTION("Pyramid quadrature values")
    {
        // Check 1 and 6 point rule
        pyramid_quadrature p1(pyramid_quadrature::point::one);
        pyramid_quadrature p8(pyramid_quadrature::point::eight);
        pyramid_quadrature p27(pyramid_quadrature::point::twenty_seven);

        REQUIRE(p1.points() == 1);
        REQUIRE(p8.points() == 8);
        REQUIRE(p27.points() == 27);

        REQUIRE(ranges::accumulate(p1.weights(), 0.0) == Approx(4.0 / 3.0));
        REQUIRE(ranges::accumulate(p8.weights(), 0.0) == Approx(4.0 / 3.0));
        REQUIRE(ranges::accumulate(p27.weights(), 0.0) == Approx(4.0 / 3.0));
    }
    SECTION("Five node evaluation")
    {
        pyramid5 element(pyramid_quadrature::point::one);

        REQUIRE(element.number_of_nodes() == 5);
        REQUIRE(element.quadrature().points() == 1);

        element.quadrature().for_each([&](auto const& femval, auto) {
            auto const& [N, rhea] = femval;

            REQUIRE(N.sum() == Approx(1.0));

            REQUIRE(rhea.col(0).sum() == Approx(0.0).margin(ZERO_MARGIN));
            REQUIRE(rhea.col(1).sum() == Approx(0.0).margin(ZERO_MARGIN));
            REQUIRE(rhea.col(2).sum() == Approx(0.0).margin(ZERO_MARGIN));
        });
        REQUIRE(element.local_quadrature_extrapolation().rows() == 5);
        REQUIRE(element.local_quadrature_extrapolation().cols() == 1);
        REQUIRE(element.local_quadrature_extrapolation().allFinite());
    }
    SECTION("Thirteen node evaluation")
    {
        pyramid13 element(pyramid_quadrature::point::eight);

        REQUIRE(element.number_of_nodes() == 13);
        REQUIRE(element.quadrature().points() == 8);

        element.quadrature().for_each([&](auto const& femval, auto) {
            auto const& [N, rhea] = femval;

            REQUIRE(N.sum() == Approx(1.0));

            REQUIRE(rhea.col(0).sum() == Approx(0.0).margin(ZERO_MARGIN));
            REQUIRE(rhea.col(1).sum() == Approx(0.0).margin(ZERO_MARGIN));
            REQUIRE(rhea.col(2).sum() == Approx(0.0).margin(ZERO_MARGIN));
        });
        REQUIRE(element.local_quadrature_extrapolation().rows() == 13);
        REQUIRE(element.local_quadrature_extrapolation().cols() == 8);
        REQUIRE(element.local_quadrature_extrapolation().allFinite());
    }
    SECTION("Five node volume - one point")
    {
        pyramid5 element(pyramid_quadrature::point::one);

        auto const base = 1.0;
        auto const height = 1.0;

        matrix3x x(3, 5);
        x << 0.0, base, base, 0.0, height / 2.0, //
            0.0, 0.0, base, base, height / 2.0,  //
            0.0, 0.0, 0.0, 0.0, height;

        REQUIRE(element.compute_measure(x) == Approx(1.0 / 3.0));
    }
    SECTION("Five node volume - eight point")
    {
        pyramid5 element(pyramid_quadrature::point::eight);

        auto constexpr base = 1.0;
        auto constexpr height = 1.0;

        matrix3x x(3, 5);
        x << 0.0, base, base, 0.0, height / 2.0, //
            0.0, 0.0, base, base, height / 2.0,  //
            0.0, 0.0, 0.0, 0.0, height;

        REQUIRE(element.compute_measure(x) == Approx(1.0 / 3.0));
    }
    SECTION("Thirteen node volume - eight point")
    {
        pyramid13 element(pyramid_quadrature::point::eight);

        matrix3x x(3, 13);
        x << 1.0, 0.0, -1.0, -1.0, -1.0, 0.0, 1.0, 1.0, 0.5, -0.5, -0.5, 0.5, 0.0, //
            1.0, 1.0, 1.0, 0.0, -1.0, -1.0, -1.0, 0.0, 0.5, 0.5, -0.5, -0.5, 0.0,  //
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.5, 0.5, 0.5, 0.5, 1.0;

        REQUIRE(element.compute_measure(x) == Approx(4.0 / 3.0));
    }
}
