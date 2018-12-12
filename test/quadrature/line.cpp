
#include <catch2/catch.hpp>

#include "quadrature/line/line_quadrature.hpp"
#include "interpolations/line.hpp"

#include "mesh/element_topology.hpp"
#include "interpolations/interpolation_factory.hpp"

#include "io/json.hpp"

#include <range/v3/numeric/accumulate.hpp>

constexpr auto ZERO_MARGIN = 1.0e-5;

using namespace neon;

json full() { return json::parse("{\"element_options\" : {\"quadrature\" : \"full\"}}"); }
json reduced() { return json::parse("{\"element_options\" : {\"quadrature\" : \"reduced\"}}"); }

TEST_CASE("Line quadrature scheme test", "[line_quadrature]")
{
    SECTION("Line Gauss Quadrature")
    {
        // Check 1 and 3 point rule
        line_quadrature l1(line_quadrature::point::one);
        line_quadrature l2(line_quadrature::point::two);
        line_quadrature l3(line_quadrature::point::three);

        REQUIRE(l1.points() == 1);
        REQUIRE(l2.points() == 2);
        REQUIRE(l3.points() == 3);

        REQUIRE(ranges::accumulate(l1.weights(), 0.0) == Approx(2.0));
        REQUIRE(ranges::accumulate(l2.weights(), 0.0) == Approx(2.0));
        REQUIRE(ranges::accumulate(l3.weights(), 0.0) == Approx(2.0));
    }
    SECTION("line 2 interpolation function - one point")
    {
        line2 line(line_quadrature::point::one);

        REQUIRE(line.number_of_nodes() == 2);
        REQUIRE(line.quadrature().points() == 1);

        line.quadrature().for_each([&](auto const& femval, auto) {
            auto const& [N, rhea] = femval;

            REQUIRE(N.sum() == Approx(1.0));

            REQUIRE(rhea.sum() == Approx(0.0).margin(ZERO_MARGIN));
        });
        REQUIRE(line.local_quadrature_extrapolation().rows() == 2);
        REQUIRE(line.local_quadrature_extrapolation().cols() == 1);
        REQUIRE(line.local_quadrature_extrapolation().allFinite());
    }
    SECTION("line 2 length")
    {
        line2 line(line_quadrature::point::one);

        matrix x(3, 2);
        x << 0.0, 1.0, //
            0.0, 0.0,  //
            0.0, 0.0;

        REQUIRE(line.compute_measure(x) == Approx(1.0));
    }
    SECTION("line 3 interpolation function - two point")
    {
        line3 line(line_quadrature::point::three);

        REQUIRE(line.number_of_nodes() == 3);
        REQUIRE(line.quadrature().points() == 3);

        line.quadrature().for_each([&](auto const& femval, auto) {
            auto const& [N, dN] = femval;

            REQUIRE(N.sum() == Approx(1.0));

            REQUIRE(dN.sum() == Approx(0.0).margin(ZERO_MARGIN));
        });
        REQUIRE(line.local_quadrature_extrapolation().rows() == 3);
        REQUIRE(line.local_quadrature_extrapolation().cols() == 3);
        REQUIRE(line.local_quadrature_extrapolation().allFinite());
    }
    SECTION("line 3 length")
    {
        line3 patch(line_quadrature::point::three);

        matrix x(3, 3);
        x << 0.0, 0.5, 1.0, //
            0.0, 0.0, 0.0,  //
            0.0, 0.0, 0.0;

        REQUIRE(patch.compute_measure(x) == Approx(1.0));
    }
    SECTION("Virtual methods check")
    {
        REQUIRE(make_line_interpolation(element_topology::line2, full())->number_of_nodes() == 2);
        REQUIRE(make_line_interpolation(element_topology::line3, full())->number_of_nodes() == 3);
    }
}
