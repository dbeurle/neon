
#include <catch2/catch.hpp>

#include "quadrature/hexahedron/hexahedron_quadrature.hpp"
#include "interpolations/hexahedron.hpp"

#include "mesh/element_topology.hpp"
#include "interpolations/interpolation_factory.hpp"

#include "io/json.hpp"

#include <range/v3/numeric/accumulate.hpp>

constexpr auto ZERO_MARGIN = 1.0e-5;

using namespace neon;

json full() { return json::parse("{\"element_options\" : {\"quadrature\" : \"full\"}}"); }
json reduced() { return json::parse("{\"element_options\" : {\"quadrature\" : \"reduced\"}}"); }

TEST_CASE("Hexahedron quadrature scheme test", "[hexahedron_quadrature]")
{
    SECTION("Hexahedron Gauss Quadrature")
    {
        // Check 1 and 8 point rule
        hexahedron_quadrature One(hexahedron_quadrature::point::one);
        hexahedron_quadrature Six(hexahedron_quadrature::point::six);
        hexahedron_quadrature Eight(hexahedron_quadrature::point::eight);
        hexahedron_quadrature TwentySeven(hexahedron_quadrature::point::twentyseven);

        REQUIRE(One.points() == 1);
        REQUIRE(Six.points() == 6);
        REQUIRE(Eight.points() == 8);
        REQUIRE(TwentySeven.points() == 27);

        REQUIRE(ranges::accumulate(One.weights(), 0.0) == Approx(8.0));
        REQUIRE(ranges::accumulate(Six.weights(), 0.0) == Approx(8.0));
        REQUIRE(ranges::accumulate(Eight.weights(), 0.0) == Approx(8.0));
        REQUIRE(ranges::accumulate(TwentySeven.weights(), 0.0) == Approx(8.0));
    }
    SECTION("hexahedron8 one Evaluation")
    {
        hexahedron8 hex8(hexahedron_quadrature::point::one);

        REQUIRE(hex8.number_of_nodes() == 8);

        // Check the quadrature is assigned correctly
        REQUIRE(hex8.quadrature().points() == 1);

        hex8.quadrature().for_each([&](auto const& femval, auto) {
            auto const& [N, rhea] = femval;

            // Check every entry
            for (int i = 0; i < 8; i++) REQUIRE(N(i) == Approx(0.125));

            REQUIRE(N.sum() == Approx(1.0));

            REQUIRE(rhea.col(0).sum() == Approx(0.0).margin(ZERO_MARGIN));
            REQUIRE(rhea.col(1).sum() == Approx(0.0).margin(ZERO_MARGIN));
            REQUIRE(rhea.col(2).sum() == Approx(0.0).margin(ZERO_MARGIN));
        });

        REQUIRE(hex8.local_quadrature_extrapolation().rows() == 8);
        REQUIRE(hex8.local_quadrature_extrapolation().cols() == 1);
        REQUIRE(hex8.local_quadrature_extrapolation().allFinite());
    }
    SECTION("hexahedron8 eight Evaluation")
    {
        hexahedron8 hex8(hexahedron_quadrature::point::eight);

        // Check the nodes are correct
        REQUIRE(hex8.number_of_nodes() == 8);

        // Check the quaduratre is assigned correctly
        REQUIRE(hex8.quadrature().points() == 8);

        hex8.quadrature().for_each([&](auto const& femval, auto) {
            auto const& [N, rhea] = femval;

            REQUIRE(N.sum() == Approx(1.0));

            REQUIRE(rhea.col(0).sum() == Approx(0.0).margin(ZERO_MARGIN));
            REQUIRE(rhea.col(1).sum() == Approx(0.0).margin(ZERO_MARGIN));
            REQUIRE(rhea.col(2).sum() == Approx(0.0).margin(ZERO_MARGIN));
        });

        REQUIRE(hex8.local_quadrature_extrapolation().rows() == 8);
        REQUIRE(hex8.local_quadrature_extrapolation().cols() == 8);
        REQUIRE(hex8.local_quadrature_extrapolation().allFinite());
    }
    SECTION("hexahedron8 volume evaluation")
    {
        hexahedron8 hex8(hexahedron_quadrature::point::eight);

        matrix x(3, 8);
        x << 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, //
            0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0,  //
            0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0;

        REQUIRE(hex8.compute_measure(x) == Approx(1.0));
    }
    SECTION("hexahedron20 one Evaluation")
    {
        hexahedron20 hex20(hexahedron_quadrature::point::one);

        // Check the nodes are correct
        REQUIRE(hex20.number_of_nodes() == 20);

        // Check the quaduratre is assigned correctly
        REQUIRE(hex20.quadrature().points() == 1);

        hex20.quadrature().for_each([&](auto const& femval, auto) {
            auto const& [N, rhea] = femval;

            REQUIRE(N.sum() == Approx(1.0));

            REQUIRE(rhea.col(0).sum() == Approx(0.0).margin(ZERO_MARGIN));
            REQUIRE(rhea.col(1).sum() == Approx(0.0).margin(ZERO_MARGIN));
            REQUIRE(rhea.col(2).sum() == Approx(0.0).margin(ZERO_MARGIN));
        });

        REQUIRE(hex20.local_quadrature_extrapolation().rows() == 20);
        REQUIRE(hex20.local_quadrature_extrapolation().cols() == 1);
        REQUIRE(hex20.local_quadrature_extrapolation().allFinite());
    }
    SECTION("hexahedron20 six Evaluation")
    {
        hexahedron20 hex20(hexahedron_quadrature::point::six);

        // Check the nodes are correct
        REQUIRE(hex20.number_of_nodes() == 20);

        // Check the quaduratre is assigned correctly
        REQUIRE(hex20.quadrature().points() == 6);

        hex20.quadrature().for_each([&](auto const& femval, auto) {
            auto const& [N, rhea] = femval;

            REQUIRE(N.sum() == Approx(1.0));

            REQUIRE(rhea.col(0).sum() == Approx(0.0).margin(ZERO_MARGIN));
            REQUIRE(rhea.col(1).sum() == Approx(0.0).margin(ZERO_MARGIN));
            REQUIRE(rhea.col(2).sum() == Approx(0.0).margin(ZERO_MARGIN));
        });

        REQUIRE(hex20.local_quadrature_extrapolation().rows() == 20);
        REQUIRE(hex20.local_quadrature_extrapolation().cols() == 6);
    }
    SECTION("hexahedron20 volume evaluation")
    {
        hexahedron20 hex20(hexahedron_quadrature::point::eight);

        // xyz coordinates of the unit cube
        matrix x(3, 20);
        x << 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.5, 1.0, 0.5, 0.0, 0.5, 1.0, 0.5, 0.0, 0.0,
            1.0, 1.0, 0.0, //
            0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, /**/ 0.0, 0.5, 1.0, 0.5, 0.0, 0.5, 1.0, 0.5,
            0.0, 0.0, 1.0, 1.0, //
            0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 0.5,
            0.5, 0.5, 0.5;

        REQUIRE(hex20.compute_measure(x) == Approx(1.0));
    }
    SECTION("hexahedron20 eight Evaluation")
    {
        hexahedron20 hex20(hexahedron_quadrature::point::eight);

        // Check the nodes are correct
        REQUIRE(hex20.number_of_nodes() == 20);

        // Check the quaduratre is assigned correctly
        REQUIRE(hex20.quadrature().points() == 8);

        hex20.quadrature().for_each([&](auto const& femval, auto) {
            auto const& [N, rhea] = femval;

            REQUIRE(N.sum() == Approx(1.0));

            REQUIRE(rhea.col(0).sum() == Approx(0.0).margin(ZERO_MARGIN));
            REQUIRE(rhea.col(1).sum() == Approx(0.0).margin(ZERO_MARGIN));
            REQUIRE(rhea.col(2).sum() == Approx(0.0).margin(ZERO_MARGIN));
        });

        REQUIRE(hex20.local_quadrature_extrapolation().rows() == 20);
        REQUIRE(hex20.local_quadrature_extrapolation().cols() == 8);
        REQUIRE(hex20.local_quadrature_extrapolation().allFinite());
    }
    SECTION("hexahedron20 TwentySevenPoint Evaluation")
    {
        hexahedron20 hex20(hexahedron_quadrature::point::twentyseven);

        // Check the nodes are correct
        REQUIRE(hex20.number_of_nodes() == 20);

        // Check the quaduratre is assigned correctly
        REQUIRE(hex20.quadrature().points() == 27);

        hex20.quadrature().for_each([&](auto const& femval, auto) {
            auto const& [N, rhea] = femval;

            REQUIRE(N.sum() == Approx(1.0));

            REQUIRE(rhea.col(0).sum() == Approx(0.0).margin(ZERO_MARGIN));
            REQUIRE(rhea.col(1).sum() == Approx(0.0).margin(ZERO_MARGIN));
            REQUIRE(rhea.col(2).sum() == Approx(0.0).margin(ZERO_MARGIN));
        });

        REQUIRE(hex20.local_quadrature_extrapolation().rows() == 20);
        REQUIRE(hex20.local_quadrature_extrapolation().cols() == 27);
        REQUIRE(hex20.local_quadrature_extrapolation().allFinite());
    }
    SECTION("hexahedron27 one Evaluation")
    {
        hexahedron27 hex27(hexahedron_quadrature::point::one);

        // Check the nodes are correct
        REQUIRE(hex27.number_of_nodes() == 27);

        // Check the quaduratre is assigned correctly
        REQUIRE(hex27.quadrature().points() == 1);

        hex27.quadrature().for_each([&](auto const& femval, auto) {
            auto const& [N, rhea] = femval;

            REQUIRE(N.sum() == Approx(1.0));

            REQUIRE(rhea.col(0).sum() == Approx(0.0).margin(ZERO_MARGIN));
            REQUIRE(rhea.col(1).sum() == Approx(0.0).margin(ZERO_MARGIN));
            REQUIRE(rhea.col(2).sum() == Approx(0.0).margin(ZERO_MARGIN));
        });

        REQUIRE(hex27.local_quadrature_extrapolation().rows() == 27);
        REQUIRE(hex27.local_quadrature_extrapolation().cols() == 1);
        REQUIRE(hex27.local_quadrature_extrapolation().allFinite());
    }
    SECTION("hexahedron27 six Evaluation")
    {
        hexahedron27 hex27(hexahedron_quadrature::point::six);

        // Check the nodes are correct
        REQUIRE(hex27.number_of_nodes() == 27);

        // Check the quaduratre is assigned correctly
        REQUIRE(hex27.quadrature().points() == 6);

        hex27.quadrature().for_each([&](auto const& femval, auto) {
            auto const& [N, rhea] = femval;

            REQUIRE(N.sum() == Approx(1.0));

            REQUIRE(rhea.col(0).sum() == Approx(0.0).margin(ZERO_MARGIN));
            REQUIRE(rhea.col(1).sum() == Approx(0.0).margin(ZERO_MARGIN));
            REQUIRE(rhea.col(2).sum() == Approx(0.0).margin(ZERO_MARGIN));
        });

        REQUIRE(hex27.local_quadrature_extrapolation().rows() == 27);
        REQUIRE(hex27.local_quadrature_extrapolation().cols() == 6);
    }
    SECTION("hexahedron27 eight Evaluation")
    {
        hexahedron27 hex27(hexahedron_quadrature::point::eight);

        // Check the nodes are correct
        REQUIRE(hex27.number_of_nodes() == 27);

        // Check the quaduratre is assigned correctly
        REQUIRE(hex27.quadrature().points() == 8);

        hex27.quadrature().for_each([&](auto const& femval, auto) {
            auto const& [N, rhea] = femval;

            REQUIRE(N.sum() == Approx(1.0));

            REQUIRE(rhea.col(0).sum() == Approx(0.0).margin(ZERO_MARGIN));
            REQUIRE(rhea.col(1).sum() == Approx(0.0).margin(ZERO_MARGIN));
            REQUIRE(rhea.col(2).sum() == Approx(0.0).margin(ZERO_MARGIN));
        });

        REQUIRE(hex27.local_quadrature_extrapolation().rows() == 27);
        REQUIRE(hex27.local_quadrature_extrapolation().cols() == 8);
        REQUIRE(hex27.local_quadrature_extrapolation().allFinite());
    }
    SECTION("hexahedron27 TwentySevenPoint Evaluation")
    {
        hexahedron27 hex27(hexahedron_quadrature::point::twentyseven);

        // Check the nodes are correct
        REQUIRE(hex27.number_of_nodes() == 27);

        // Check the quaduratre is assigned correctly
        REQUIRE(hex27.quadrature().points() == 27);

        hex27.quadrature().for_each([&](auto const& femval, auto) {
            auto const& [N, rhea] = femval;

            REQUIRE(N.sum() == Approx(1.0));

            REQUIRE(rhea.col(0).sum() == Approx(0.0).margin(ZERO_MARGIN));
            REQUIRE(rhea.col(1).sum() == Approx(0.0).margin(ZERO_MARGIN));
            REQUIRE(rhea.col(2).sum() == Approx(0.0).margin(ZERO_MARGIN));
        });

        REQUIRE(hex27.local_quadrature_extrapolation().rows() == 27);
        REQUIRE(hex27.local_quadrature_extrapolation().cols() == 27);
    }
    SECTION("Virtual methods check")
    {
        REQUIRE(make_volume_interpolation(element_topology::hexahedron8, full())->number_of_nodes()
                == 8);
        REQUIRE(make_volume_interpolation(element_topology::hexahedron20, full())->number_of_nodes()
                == 20);
        REQUIRE(make_volume_interpolation(element_topology::hexahedron27, full())->number_of_nodes()
                == 27);
    }
}
TEST_CASE("Hexahedron volume evaluation")
{
    SECTION("hexahedron27 volume evaluation")
    {
        SECTION("Six point rule")
        {
            hexahedron27 hex27(hexahedron_quadrature::point::six);

            // xyz coordinates of the unit cube
            matrix x(3, 27);
            x << 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, /**/ 0.5, 1.0, 0.5, 0.0, 0.5, 1.0, 0.5,
                0.0, 0.0, 1.0, 1.0, 0.0, /**/ 0.5, 0.5, 0.5, 0.5, 0.0, 1.0, 0.5, //
                0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, /**/ 0.0, 0.5, 1.0, 0.5, 0.0, 0.5, 1.0, 0.5,
                0.0, 0.0, 1.0, 1.0, /**/ 0.5, 0.5, 0.0, 1.0, 0.5, 0.5, 0.5, //
                0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, /**/ 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0,
                0.5, 0.5, 0.5, 0.5, /**/ 0.0, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5;

            REQUIRE(hex27.compute_measure(x) == Approx(1.0));

            x *= 2.0;

            REQUIRE(hex27.compute_measure(x) == Approx(8.0));
        }
        SECTION("Eight point rule")
        {
            hexahedron27 hex27(hexahedron_quadrature::point::eight);

            // xyz coordinates of the unit cube
            matrix x(3, 27);
            x << 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, /**/ 0.5, 1.0, 0.5, 0.0, 0.5, 1.0, 0.5,
                0.0, 0.0, 1.0, 1.0, 0.0, /**/ 0.5, 0.5, 0.5, 0.5, 0.0, 1.0, 0.5, //
                0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, /**/ 0.0, 0.5, 1.0, 0.5, 0.0, 0.5, 1.0, 0.5,
                0.0, 0.0, 1.0, 1.0, /**/ 0.5, 0.5, 0.0, 1.0, 0.5, 0.5, 0.5, //
                0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, /**/ 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0,
                0.5, 0.5, 0.5, 0.5, /**/ 0.0, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5;

            REQUIRE(hex27.compute_measure(x) == Approx(1.0));

            x *= 2.0;

            REQUIRE(hex27.compute_measure(x) == Approx(8.0));
        }
    }
}
