
#include <catch2/catch.hpp>

#include "quadrature/prism/prism_quadrature.hpp"
#include "interpolations/prism.hpp"

#include "interpolations/interpolation_factory.hpp"
#include "mesh/element_topology.hpp"

#include "io/json.hpp"

#include <range/v3/numeric/accumulate.hpp>

using namespace neon;

json full() { return json::parse("{\"element_options\" : {\"quadrature\" : \"full\"}}"); }
json reduced() { return json::parse("{\"element_options\" : {\"quadrature\" : \"reduced\"}}"); }

constexpr auto ZERO_MARGIN = 1.0e-5;

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

        matrix3x x(3, 6);
        x << 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, //
            0.0, 1.0, 0.0, 0.0, 1.0, 0.0,  //
            -1.0, -1.0, -1.0, 1.0, 1.0, 1.0;

        REQUIRE(pri6.compute_measure(x) == Approx(1.0));
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

        matrix3x x(3, 15);
        x << 1.0, 0.0, 0.0, 0.5, 0.0, 0.5, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.5, 0.0, 0.5, //
            0.0, 1.0, 0.0, 0.5, 0.5, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.5, 0.5, 0.0,  //
            -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0;

        REQUIRE(pri15.compute_measure(x) == Approx(1.0));
    }
    SECTION("Virtual methods check")
    {
        REQUIRE(make_volume_interpolation(element_topology::prism6, full())->number_of_nodes() == 6);
        REQUIRE(make_volume_interpolation(element_topology::prism15, full())->number_of_nodes() == 15);
    }
}
