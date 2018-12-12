
#include <catch2/catch.hpp>

#include "quadrature/quadrilateral/quadrilateral_quadrature.hpp"
#include "interpolations/quadrilateral.hpp"

#include "mesh/element_topology.hpp"
#include "interpolations/interpolation_factory.hpp"

#include "io/json.hpp"

#include <range/v3/numeric/accumulate.hpp>

constexpr auto ZERO_MARGIN = 1.0e-5;

using namespace neon;

json full() { return json::parse("{\"element_options\" : {\"quadrature\" : \"full\"}}"); }
json reduced() { return json::parse("{\"element_options\" : {\"quadrature\" : \"reduced\"}}"); }

TEST_CASE("Quadrilateral quadrature scheme test", "[quadrilateral_quadrature]")
{
    SECTION("Quadrilateral Gauss Quadrature")
    {
        // Check 1 and 8 point rule
        quadrilateral_quadrature q1(quadrilateral_quadrature::point::one);
        quadrilateral_quadrature q4(quadrilateral_quadrature::point::four);
        quadrilateral_quadrature q9(quadrilateral_quadrature::point::nine);

        REQUIRE(q1.points() == 1);
        REQUIRE(q4.points() == 4);
        REQUIRE(q9.points() == 9);

        REQUIRE(ranges::accumulate(q1.weights(), 0.0) == Approx(4.0));
        REQUIRE(ranges::accumulate(q4.weights(), 0.0) == Approx(4.0));
        REQUIRE(ranges::accumulate(q9.weights(), 0.0) == Approx(4.0));
    }
    SECTION("quadrilateral4 interpolation function - one point")
    {
        quadrilateral4 quad4(quadrilateral_quadrature::point::one);

        REQUIRE(quad4.number_of_nodes() == 4);
        REQUIRE(quad4.quadrature().points() == 1);
        REQUIRE(quad4.local_quadrature_extrapolation().rows() == 4);
        REQUIRE(quad4.local_quadrature_extrapolation().cols() == 1);

        quad4.quadrature().for_each([&](auto const& femval, auto) {
            auto const& [N, rhea] = femval;

            REQUIRE(N.sum() == Approx(1.0));

            REQUIRE(rhea.col(0).sum() == Approx(0.0).margin(ZERO_MARGIN));
            REQUIRE(rhea.col(1).sum() == Approx(0.0).margin(ZERO_MARGIN));
        });

        REQUIRE(quad4.local_quadrature_extrapolation().allFinite());
    }
    SECTION("quadrilateral4 surface area - one point")
    {
        quadrilateral4 quad4(quadrilateral_quadrature::point::one);

        matrix x(3, 4);
        x << 0.0, 1.0, 1.0, 0.0, //
            0.0, 0.0, 1.0, 1.0,  //
            0.0, 0.0, 0.0, 0.0;

        REQUIRE(quad4.compute_measure(x) == Approx(1.0));
    }
    SECTION("quadrilateral4 interpolation function - four point")
    {
        quadrilateral4 quad4(quadrilateral_quadrature::point::four);

        REQUIRE(quad4.number_of_nodes() == 4);
        REQUIRE(quad4.quadrature().points() == 4);
        REQUIRE(quad4.local_quadrature_extrapolation().rows() == 4);
        REQUIRE(quad4.local_quadrature_extrapolation().cols() == 4);

        quad4.quadrature().for_each([&](auto const& femval, auto) {
            auto const& [N, rhea] = femval;

            REQUIRE(N.sum() == Approx(1.0));

            REQUIRE(rhea.col(0).sum() == Approx(0.0).margin(ZERO_MARGIN));
            REQUIRE(rhea.col(1).sum() == Approx(0.0).margin(ZERO_MARGIN));
        });

        REQUIRE(quad4.local_quadrature_extrapolation().allFinite());
    }
    SECTION("quadrilateral4 surface area - four point")
    {
        quadrilateral4 quad4(quadrilateral_quadrature::point::four);

        matrix x(3, 4);
        x << 0.0, 1.0, 1.0, 0.0, //
            0.0, 0.0, 1.0, 1.0,  //
            0.0, 0.0, 0.0, 0.0;

        REQUIRE(quad4.compute_measure(x) == Approx(1.0));
    }
    SECTION("quadrilateral8 interpolation function - four point")
    {
        quadrilateral8 quad8(quadrilateral_quadrature::point::four);

        REQUIRE(quad8.number_of_nodes() == 8);
        REQUIRE(quad8.quadrature().points() == 4);
        REQUIRE(quad8.local_quadrature_extrapolation().rows() == 8);
        REQUIRE(quad8.local_quadrature_extrapolation().cols() == 4);

        quad8.quadrature().for_each([&](auto const& femval, auto) {
            auto const& [N, rhea] = femval;

            REQUIRE(N.sum() == Approx(1.0));

            REQUIRE(rhea.col(0).sum() == Approx(0.0).margin(ZERO_MARGIN));
            REQUIRE(rhea.col(1).sum() == Approx(0.0).margin(ZERO_MARGIN));
        });

        REQUIRE(quad8.local_quadrature_extrapolation().allFinite());
    }
    SECTION("quadrilateral8 surface area - four point")
    {
        quadrilateral8 quad8(quadrilateral_quadrature::point::four);

        matrix x(3, 8);
        x << 0.0, 1.0, 1.0, 0.0, 0.5, 1.0, 0.5, 0.0, //
            0.0, 0.0, 1.0, 1.0, 0.0, 0.5, 1.0, 0.5,  //
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0;

        REQUIRE(quad8.compute_measure(x) == Approx(1.0));
    }
    SECTION("quadrilateral8 interpolation function - nine point")
    {
        quadrilateral8 quad8(quadrilateral_quadrature::point::nine);

        REQUIRE(quad8.number_of_nodes() == 8);
        REQUIRE(quad8.quadrature().points() == 9);
        REQUIRE(quad8.local_quadrature_extrapolation().rows() == 8);
        REQUIRE(quad8.local_quadrature_extrapolation().cols() == 9);

        quad8.quadrature().for_each([&](auto const& femval, auto) {
            auto const& [N, rhea] = femval;

            REQUIRE(N.sum() == Approx(1.0));

            REQUIRE(rhea.col(0).sum() == Approx(0.0).margin(ZERO_MARGIN));
            REQUIRE(rhea.col(1).sum() == Approx(0.0).margin(ZERO_MARGIN));
        });
        REQUIRE(quad8.local_quadrature_extrapolation().allFinite());
    }
    SECTION("quadrilateral8 interpolation function - nine point")
    {
        quadrilateral8 quad8(quadrilateral_quadrature::point::nine);

        matrix x(3, 8);
        x << 0.0, 1.0, 1.0, 0.0, 0.5, 1.0, 0.5, 0.0, //
            0.0, 0.0, 1.0, 1.0, 0.0, 0.5, 1.0, 0.5,  //
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0;

        REQUIRE(quad8.compute_measure(x) == Approx(1.0));
    }
    SECTION("quadrilateral9 interpolation function - four point")
    {
        quadrilateral9 quad9(quadrilateral_quadrature::point::four);

        REQUIRE(quad9.number_of_nodes() == 9);
        REQUIRE(quad9.quadrature().points() == 4);
        REQUIRE(quad9.local_quadrature_extrapolation().rows() == 9);
        REQUIRE(quad9.local_quadrature_extrapolation().cols() == 4);

        quad9.quadrature().for_each([&](auto const& femval, auto) {
            auto const& [N, rhea] = femval;

            REQUIRE(N.sum() == Approx(1.0));

            REQUIRE(rhea.col(0).sum() == Approx(0.0).margin(ZERO_MARGIN));
            REQUIRE(rhea.col(1).sum() == Approx(0.0).margin(ZERO_MARGIN));
        });
        REQUIRE(quad9.local_quadrature_extrapolation().allFinite());
    }
    SECTION("quadrilateral9 surface area - four point")
    {
        quadrilateral9 quad9(quadrilateral_quadrature::point::four);

        matrix x(3, 9);
        x << 0.0, 1.0, 1.0, 0.0, 0.5, 1.0, 0.5, 0.0, 0.5, //
            0.0, 0.0, 1.0, 1.0, 0.0, 0.5, 1.0, 0.5, 0.5,  //
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0;

        REQUIRE(quad9.compute_measure(x) == Approx(1.0));
    }
    SECTION("quadrilateral9 interpolation function - nine point")
    {
        quadrilateral9 quad9(quadrilateral_quadrature::point::nine);

        REQUIRE(quad9.number_of_nodes() == 9);
        REQUIRE(quad9.quadrature().points() == 9);
        REQUIRE(quad9.local_quadrature_extrapolation().rows() == 9);
        REQUIRE(quad9.local_quadrature_extrapolation().cols() == 9);

        quad9.quadrature().for_each([&](auto const& femval, auto) {
            auto const& [N, rhea] = femval;

            REQUIRE(N.sum() == Approx(1.0));

            REQUIRE(rhea.col(0).sum() == Approx(0.0).margin(ZERO_MARGIN));
            REQUIRE(rhea.col(1).sum() == Approx(0.0).margin(ZERO_MARGIN));
        });
        REQUIRE(quad9.local_quadrature_extrapolation().allFinite());
    }
    SECTION("quadrilateral9 interpolation function - nine point")
    {
        quadrilateral9 quad9(quadrilateral_quadrature::point::nine);

        matrix x(3, 9);
        x << 0.0, 1.0, 1.0, 0.0, 0.5, 1.0, 0.5, 0.0, 0.5, //
            0.0, 0.0, 1.0, 1.0, 0.0, 0.5, 1.0, 0.5, 0.5,  //
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0;

        REQUIRE(quad9.compute_measure(x) == Approx(1.0));
    }
    SECTION("Virtual methods check")
    {
        REQUIRE(make_surface_interpolation(element_topology::quadrilateral4, full())->number_of_nodes()
                == 4);
        REQUIRE(make_surface_interpolation(element_topology::quadrilateral8, full())->number_of_nodes()
                == 8);
        REQUIRE(make_surface_interpolation(element_topology::quadrilateral9, full())->number_of_nodes()
                == 9);
    }
}
