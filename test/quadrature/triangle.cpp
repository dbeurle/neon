
#include <catch2/catch.hpp>

#include "quadrature/triangle/triangle_quadrature.hpp"
#include "interpolations/triangle.hpp"

#include "mesh/element_topology.hpp"
#include "interpolations/interpolation_factory.hpp"

#include "io/json.hpp"

#include <range/v3/numeric/accumulate.hpp>

constexpr auto ZERO_MARGIN = 1.0e-5;

using namespace neon;

json full() { return json::parse("{\"element_options\" : {\"quadrature\" : \"full\"}}"); }
json reduced() { return json::parse("{\"element_options\" : {\"quadrature\" : \"reduced\"}}"); }

TEST_CASE("Triangle quadrature scheme test", "[triangle_quadrature]")
{
    SECTION("triangle Gauss quadrature")
    {
        // Check 1 and 8 point rule
        triangle_quadrature t1(triangle_quadrature::point::one);
        triangle_quadrature t3(triangle_quadrature::point::three);
        triangle_quadrature t4(triangle_quadrature::point::four);

        REQUIRE(t1.points() == 1);
        REQUIRE(t3.points() == 3);
        REQUIRE(t4.points() == 4);

        REQUIRE(ranges::accumulate(t1.weights(), 0.0) == Approx(0.5));
        REQUIRE(ranges::accumulate(t3.weights(), 0.0) == Approx(0.5));
        REQUIRE(ranges::accumulate(t4.weights(), 0.0) == Approx(0.5));
    }
    SECTION("triangle3 interpolation function - one point")
    {
        triangle3 tri3(triangle_quadrature::point::one);

        REQUIRE(tri3.number_of_nodes() == 3);
        REQUIRE(tri3.quadrature().points() == 1);
        REQUIRE(tri3.local_quadrature_extrapolation().rows() == 3);
        REQUIRE(tri3.local_quadrature_extrapolation().cols() == 1);

        tri3.quadrature().for_each([&](auto const& femval, auto) {
            auto const& [N, rhea] = femval;

            REQUIRE(N.sum() == Approx(1.0));

            REQUIRE(rhea.col(0).sum() == Approx(0.0).margin(ZERO_MARGIN));
            REQUIRE(rhea.col(1).sum() == Approx(0.0).margin(ZERO_MARGIN));
        });
        REQUIRE(tri3.local_quadrature_extrapolation().allFinite());
    }
    SECTION("triangle6 interpolation function - three point")
    {
        triangle6 tri6(triangle_quadrature::point::three);

        REQUIRE(tri6.number_of_nodes() == 6);
        REQUIRE(tri6.quadrature().points() == 3);

        tri6.quadrature().for_each([&](auto const& femval, auto) {
            auto const& [N, rhea] = femval;

            REQUIRE(N.sum() == Approx(1.0));

            REQUIRE(rhea.col(0).sum() == Approx(0.0).margin(ZERO_MARGIN));
            REQUIRE(rhea.col(1).sum() == Approx(0.0).margin(ZERO_MARGIN));
        });
        REQUIRE(tri6.local_quadrature_extrapolation().allFinite());
    }
    SECTION("triangle6 interpolation function - four point")
    {
        triangle6 tri6(triangle_quadrature::point::four);

        REQUIRE(tri6.number_of_nodes() == 6);
        REQUIRE(tri6.quadrature().points() == 4);

        tri6.quadrature().for_each([&](auto const& femval, auto) {
            auto const& [N, rhea] = femval;

            REQUIRE(N.sum() == Approx(1.0));

            REQUIRE(rhea.col(0).sum() == Approx(0.0).margin(ZERO_MARGIN));
            REQUIRE(rhea.col(1).sum() == Approx(0.0).margin(ZERO_MARGIN));
        });
        REQUIRE(tri6.local_quadrature_extrapolation().allFinite());
    }
    SECTION("triangle3 surface area")
    {
        triangle3 tri3(triangle_quadrature::point::one);

        matrix x(3, 3);
        x << 0.0, 0.0, 0.0, //
            1.0, 0.0, 0.0,  //
            0.0, 1.0, 0.0;
        x.transposeInPlace();

        REQUIRE(tri3.compute_measure(x) == Approx(0.5));
    }
    SECTION("triangle6 surface area")
    {
        triangle6 patch(triangle_quadrature::point::three);

        matrix x(6, 3);
        x << 0.0, 0.0, 0.0, //
            1.0, 0.0, 0.0,  //
            0.0, 1.0, 0.0,  //
            0.5, 0.0, 0.0,  //
            0.5, 0.5, 0.0,  //
            0.0, 0.5, 0.0;
        x.transposeInPlace();

        REQUIRE(patch.compute_measure(x) == Approx(0.5));
    }
    SECTION("Virtual methods check")
    {
        REQUIRE(make_surface_interpolation(element_topology::triangle3, full())->number_of_nodes()
                == 3);
        REQUIRE(make_surface_interpolation(element_topology::triangle6, full())->number_of_nodes()
                == 6);
    }
}
