
#include <catch.hpp>

#include "quadrature/hexahedron_quadrature.hpp"
#include "quadrature/line_quadrature.hpp"
#include "quadrature/prism_quadrature.hpp"
#include "quadrature/quadrilateral_quadrature.hpp"
#include "quadrature/tetrahedron_quadrature.hpp"
#include "quadrature/triangle_quadrature.hpp"
#include "quadrature/unit_sphere_quadrature.hpp"

#include "interpolations/line.hpp"
#include "interpolations/quadrilateral.hpp"
#include "interpolations/hexahedron.hpp"
#include "interpolations/tetrahedron.hpp"
#include "interpolations/triangle.hpp"
#include "interpolations/prism.hpp"

#include "interpolations/interpolation_factory.hpp"
#include "mesh/element_topology.hpp"

#include "io/json.hpp"

#include <range/v3/numeric/accumulate.hpp>

using namespace neon;

json full() { return json::parse("{\"ElementOptions\" : {\"Quadrature\" : \"Full\"}}"); }
json reduced() { return json::parse("{\"ElementOptions\" : {\"Quadrature\" : \"Reduced\"}}"); }

constexpr auto ZERO_MARGIN = 1.0e-5;

TEST_CASE("Line quadrature scheme test", "[line_quadrature]")
{
    SECTION("Line Gauss Quadrature")
    {
        // Check 1 and 8 point rule
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

        REQUIRE(line.nodes() == 2);
        REQUIRE(line.quadrature().points() == 1);

        line.quadrature().for_each([&](auto const& femval, auto const& l) {
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
        line2 patch(line_quadrature::point::one);

        matrix x(3, 2);
        x << 0.0, 1.0, //
            0.0, 0.0,  //
            0.0, 0.0;

        REQUIRE(patch.compute_measure(x) == Approx(1.0));
    }
    SECTION("line 3 interpolation function - two point")
    {
        line3 line(line_quadrature::point::two);

        REQUIRE(line.nodes() == 3);
        REQUIRE(line.quadrature().points() == 2);

        line.quadrature().for_each([&](auto const& femval, auto const& l) {
            auto const& [N, rhea] = femval;

            REQUIRE(N.sum() == Approx(1.0));

            REQUIRE(rhea.sum() == Approx(0.0).margin(ZERO_MARGIN));
        });
        REQUIRE(line.local_quadrature_extrapolation().rows() == 3);
        REQUIRE(line.local_quadrature_extrapolation().cols() == 2);
        REQUIRE(line.local_quadrature_extrapolation().allFinite());
    }
    SECTION("line 3 length")
    {
        line3 patch(line_quadrature::point::two);

        matrix x(3, 3);
        x << 0.0, 0.5, 1.0, //
            0.0, 0.0, 0.0,  //
            0.0, 0.0, 0.0;

        REQUIRE(patch.compute_measure(x) == Approx(1.0));
    }
    SECTION("Virtual methods check")
    {
        REQUIRE(make_line_interpolation(element_topology::line2, full())->nodes() == 2);
        REQUIRE(make_line_interpolation(element_topology::line3, full())->nodes() == 3);
    }
}
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

        REQUIRE(quad4.nodes() == 4);
        REQUIRE(quad4.quadrature().points() == 1);
        REQUIRE(quad4.local_quadrature_extrapolation().rows() == 4);
        REQUIRE(quad4.local_quadrature_extrapolation().cols() == 1);

        quad4.quadrature().for_each([&](auto const& femval, auto const& l) {
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

        REQUIRE(quad4.nodes() == 4);
        REQUIRE(quad4.quadrature().points() == 4);
        REQUIRE(quad4.local_quadrature_extrapolation().rows() == 4);
        REQUIRE(quad4.local_quadrature_extrapolation().cols() == 4);

        quad4.quadrature().for_each([&](auto const& femval, auto const& l) {
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

        REQUIRE(quad8.nodes() == 8);
        REQUIRE(quad8.quadrature().points() == 4);
        REQUIRE(quad8.local_quadrature_extrapolation().rows() == 8);
        REQUIRE(quad8.local_quadrature_extrapolation().cols() == 4);

        quad8.quadrature().for_each([&](auto const& femval, auto const& l) {
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

        REQUIRE(quad8.nodes() == 8);
        REQUIRE(quad8.quadrature().points() == 9);
        REQUIRE(quad8.local_quadrature_extrapolation().rows() == 8);
        REQUIRE(quad8.local_quadrature_extrapolation().cols() == 9);

        quad8.quadrature().for_each([&](auto const& femval, auto const& l) {
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

        REQUIRE(quad9.nodes() == 9);
        REQUIRE(quad9.quadrature().points() == 4);
        REQUIRE(quad9.local_quadrature_extrapolation().rows() == 9);
        REQUIRE(quad9.local_quadrature_extrapolation().cols() == 4);

        quad9.quadrature().for_each([&](auto const& femval, auto const& l) {
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

        REQUIRE(quad9.nodes() == 9);
        REQUIRE(quad9.quadrature().points() == 9);
        REQUIRE(quad9.local_quadrature_extrapolation().rows() == 9);
        REQUIRE(quad9.local_quadrature_extrapolation().cols() == 9);

        quad9.quadrature().for_each([&](auto const& femval, auto const& l) {
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
        REQUIRE(make_surface_interpolation(element_topology::quadrilateral4, full())->nodes() == 4);
        REQUIRE(make_surface_interpolation(element_topology::quadrilateral8, full())->nodes() == 8);
        REQUIRE(make_surface_interpolation(element_topology::quadrilateral9, full())->nodes() == 9);
    }
}
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

        REQUIRE(tri3.nodes() == 3);
        REQUIRE(tri3.quadrature().points() == 1);
        REQUIRE(tri3.local_quadrature_extrapolation().rows() == 3);
        REQUIRE(tri3.local_quadrature_extrapolation().cols() == 1);

        tri3.quadrature().for_each([&](auto const& femval, auto const& l) {
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

        REQUIRE(tri6.nodes() == 6);
        REQUIRE(tri6.quadrature().points() == 3);

        tri6.quadrature().for_each([&](auto const& femval, auto const& l) {
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

        REQUIRE(tri6.nodes() == 6);
        REQUIRE(tri6.quadrature().points() == 4);

        tri6.quadrature().for_each([&](auto const& femval, auto const& l) {
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
        REQUIRE(make_surface_interpolation(element_topology::triangle3, full())->nodes() == 3);
        REQUIRE(make_surface_interpolation(element_topology::triangle6, full())->nodes() == 6);
    }
}
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

        REQUIRE(hex8.nodes() == 8);

        // Check the quadrature is assigned correctly
        REQUIRE(hex8.quadrature().points() == 1);

        hex8.quadrature().for_each([&](auto const& femval, auto const& l) {
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
        REQUIRE(hex8.nodes() == 8);

        // Check the quaduratre is assigned correctly
        REQUIRE(hex8.quadrature().points() == 8);

        hex8.quadrature().for_each([&](auto const& femval, auto const& l) {
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
        REQUIRE(hex20.nodes() == 20);

        // Check the quaduratre is assigned correctly
        REQUIRE(hex20.quadrature().points() == 1);

        hex20.quadrature().for_each([&](auto const& femval, auto const& l) {
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
        REQUIRE(hex20.nodes() == 20);

        // Check the quaduratre is assigned correctly
        REQUIRE(hex20.quadrature().points() == 6);

        hex20.quadrature().for_each([&](auto const& femval, auto const& l) {
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
        REQUIRE(hex20.nodes() == 20);

        // Check the quaduratre is assigned correctly
        REQUIRE(hex20.quadrature().points() == 8);

        hex20.quadrature().for_each([&](auto const& femval, auto const& l) {
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
        REQUIRE(hex20.nodes() == 20);

        // Check the quaduratre is assigned correctly
        REQUIRE(hex20.quadrature().points() == 27);

        hex20.quadrature().for_each([&](auto const& femval, auto const& l) {
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
        REQUIRE(hex27.nodes() == 27);

        // Check the quaduratre is assigned correctly
        REQUIRE(hex27.quadrature().points() == 1);

        hex27.quadrature().for_each([&](auto const& femval, auto const& l) {
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
        REQUIRE(hex27.nodes() == 27);

        // Check the quaduratre is assigned correctly
        REQUIRE(hex27.quadrature().points() == 6);

        hex27.quadrature().for_each([&](auto const& femval, auto const& l) {
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
        REQUIRE(hex27.nodes() == 27);

        // Check the quaduratre is assigned correctly
        REQUIRE(hex27.quadrature().points() == 8);

        hex27.quadrature().for_each([&](auto const& femval, auto const& l) {
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
        REQUIRE(hex27.nodes() == 27);

        // Check the quaduratre is assigned correctly
        REQUIRE(hex27.quadrature().points() == 27);

        hex27.quadrature().for_each([&](auto const& femval, auto const& l) {
            auto const& [N, rhea] = femval;

            REQUIRE(N.sum() == Approx(1.0));

            REQUIRE(rhea.col(0).sum() == Approx(0.0).margin(ZERO_MARGIN));
            REQUIRE(rhea.col(1).sum() == Approx(0.0).margin(ZERO_MARGIN));
            REQUIRE(rhea.col(2).sum() == Approx(0.0).margin(ZERO_MARGIN));
        });

        REQUIRE(hex27.local_quadrature_extrapolation().rows() == 27);
        REQUIRE(hex27.local_quadrature_extrapolation().cols() == 27);
    }
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
    SECTION("Virtual methods check")
    {
        REQUIRE(make_volume_interpolation(element_topology::hexahedron8, full())->nodes() == 8);
        REQUIRE(make_volume_interpolation(element_topology::hexahedron20, full())->nodes() == 20);
        REQUIRE(make_volume_interpolation(element_topology::hexahedron27, full())->nodes() == 27);
    }
}
TEST_CASE("Tetrahedron quadrature scheme test", "[tetrahedron_quadrature]")
{
    SECTION("Tetrahedron Gauss Quadrature")
    {
        tetrahedron_quadrature One(tetrahedron_quadrature::point::one);
        tetrahedron_quadrature Four(tetrahedron_quadrature::point::four);
        tetrahedron_quadrature Five(tetrahedron_quadrature::point::five);

        // Check the number of quadrature points are correctly set
        REQUIRE(One.points() == 1);
        REQUIRE(Four.points() == 4);
        REQUIRE(Five.points() == 5);

        // Check the weightings add to 1/6 for each scheme
        REQUIRE(ranges::accumulate(One.weights(), 0.0) == Approx(1.0 / 6.0));
        REQUIRE(ranges::accumulate(Four.weights(), 0.0) == Approx(1.0 / 6.0));
        REQUIRE(ranges::accumulate(Five.weights(), 0.0) == Approx(1.0 / 6.0));
    }
    SECTION("tetrahedron4 one Evaluation")
    {
        tetrahedron4 tet4(tetrahedron_quadrature::point::one);

        REQUIRE(tet4.nodes() == 4);
        REQUIRE(tet4.quadrature().points() == 1);

        tet4.quadrature().for_each([&](auto const& femval, auto const& l) {
            auto const& [N, rhea] = femval;

            REQUIRE(N.sum() == Approx(1.0));

            REQUIRE(rhea.col(0).sum() == Approx(0.0).margin(ZERO_MARGIN));
            REQUIRE(rhea.col(1).sum() == Approx(0.0).margin(ZERO_MARGIN));
            REQUIRE(rhea.col(2).sum() == Approx(0.0).margin(ZERO_MARGIN));
        });

        REQUIRE(tet4.local_quadrature_extrapolation().rows() == 4);
        REQUIRE(tet4.local_quadrature_extrapolation().cols() == 1);
        REQUIRE(tet4.local_quadrature_extrapolation().allFinite());
    }
    SECTION("tetrahedron10 one Evaluation")
    {
        tetrahedron10 tet10(tetrahedron_quadrature::point::one);

        REQUIRE(tet10.nodes() == 10);

        REQUIRE(tet10.quadrature().points() == 1);

        tet10.quadrature().for_each([&](auto const& femval, auto const& l) {
            auto const& [N, rhea] = femval;

            REQUIRE(N.sum() == Approx(1.0));

            REQUIRE(rhea.col(0).sum() == Approx(0.0).margin(ZERO_MARGIN));
            REQUIRE(rhea.col(1).sum() == Approx(0.0).margin(ZERO_MARGIN));
            REQUIRE(rhea.col(2).sum() == Approx(0.0).margin(ZERO_MARGIN));
        });

        REQUIRE(tet10.local_quadrature_extrapolation().rows() == 10);
        REQUIRE(tet10.local_quadrature_extrapolation().cols() == 1);
        REQUIRE(tet10.local_quadrature_extrapolation().allFinite());
    }
    SECTION("tetrahedron10 four Evaluation")
    {
        tetrahedron10 tet10(tetrahedron_quadrature::point::four);

        REQUIRE(tet10.nodes() == 10);

        REQUIRE(tet10.quadrature().points() == 4);

        tet10.quadrature().for_each([&](auto const& femval, auto const& l) {
            auto const& [N, rhea] = femval;

            REQUIRE(N.sum() == Approx(1.0));

            REQUIRE(rhea.col(0).sum() == Approx(0.0).margin(ZERO_MARGIN));
            REQUIRE(rhea.col(1).sum() == Approx(0.0).margin(ZERO_MARGIN));
            REQUIRE(rhea.col(2).sum() == Approx(0.0).margin(ZERO_MARGIN));
        });

        REQUIRE(tet10.local_quadrature_extrapolation().rows() == 10);
        REQUIRE(tet10.local_quadrature_extrapolation().cols() == 4);
        REQUIRE(tet10.local_quadrature_extrapolation().allFinite());
    }
    SECTION("tetrahedron10 five Evaluation")
    {
        tetrahedron10 tet10(tetrahedron_quadrature::point::five);

        REQUIRE(tet10.nodes() == 10);

        REQUIRE(tet10.quadrature().points() == 5);

        tet10.quadrature().for_each([&](auto const& femval, auto const& l) {
            auto const& [N, rhea] = femval;

            REQUIRE(N.sum() == Approx(1.0));

            REQUIRE(rhea.col(0).sum() == Approx(0.0).margin(ZERO_MARGIN));
            REQUIRE(rhea.col(1).sum() == Approx(0.0).margin(ZERO_MARGIN));
            REQUIRE(rhea.col(2).sum() == Approx(0.0).margin(ZERO_MARGIN));
        });

        REQUIRE(tet10.local_quadrature_extrapolation().rows() == 10);
        REQUIRE(tet10.local_quadrature_extrapolation().cols() == 5);
        REQUIRE(tet10.local_quadrature_extrapolation().allFinite());
    }
    SECTION("tetrahedron10 Jacobian and volume check")
    {
        SECTION("Parent element")
        {
            // Check that we have a Jacobian of one with a parent element
            matrix x(10, 3);
            x << 1.0, 0.0, 0.0, //
                0.0, 1.0, 0.0,  //
                0.0, 0.0, 1.0,  //
                0.0, 0.0, 0.0,  //
                0.5, 0.5, 0.0,  //
                0.0, 0.5, 0.5,  //
                0.0, 0.0, 0.5,  //
                0.5, 0.0, 0.0,  //
                0.5, 0.0, 0.5,  //
                0.0, 0.5, 0.0;
            x.transposeInPlace();

            SECTION("One point")
            {
                tetrahedron10 tet10(tetrahedron_quadrature::point::one);

                tet10.quadrature().for_each([&](auto const& femval, auto const& l) {
                    auto const& [N, rhea] = femval;

                    matrix3 const Jacobian = x * rhea;

                    REQUIRE(Jacobian.determinant() == Approx(1.0));
                });
                double vol = tet10.quadrature().integrate(0.0, [&](auto const& femval, auto const& l) {
                    auto const& [N, rhea] = femval;

                    matrix3 const Jacobian = x * rhea;

                    return Jacobian.determinant();
                });
                REQUIRE(vol == Approx(1.0 / 6.0));
            }
            SECTION("Four point")
            {
                tetrahedron10 tet10(tetrahedron_quadrature::point::four);

                tet10.quadrature().for_each([&](auto const& femval, auto const& l) {
                    auto const& [N, rhea] = femval;

                    matrix3 const Jacobian = x * rhea;

                    REQUIRE(Jacobian.determinant() == Approx(1.0));
                });

                double vol = tet10.quadrature().integrate(0.0, [&](auto const& femval, auto const& l) {
                    auto const& [N, rhea] = femval;

                    matrix3 const Jacobian = x * rhea;

                    return Jacobian.determinant();
                });
                REQUIRE(vol == Approx(1.0 / 6.0));
            }
            SECTION("Five point")
            {
                tetrahedron10 tet10(tetrahedron_quadrature::point::five);

                tet10.quadrature().for_each([&](auto const& femval, auto const& l) {
                    auto const& [N, rhea] = femval;

                    matrix3 const Jacobian = x * rhea;

                    REQUIRE(Jacobian.determinant() == Approx(1.0));
                });

                double vol = tet10.quadrature().integrate(0.0, [&](auto const& femval, auto const& l) {
                    auto const& [N, rhea] = femval;

                    matrix3 const Jacobian = x * rhea;

                    return Jacobian.determinant();
                });
                REQUIRE(vol == Approx(1.0 / 6.0));
            }
        }
        SECTION("Volume of tetrahedron with unit length sides")
        {
            tetrahedron10 tet10(tetrahedron_quadrature::point::four);

            matrix x_vertices(4, 3);
            x_vertices << 1.0, 0.0, 0.0,                          // Node 0
                0.5, 0.5 * std::sqrt(3.0), 0.0,                   // Node 1
                0.5, std::sqrt(1.0 / 12.0), std::sqrt(2.0 / 3.0), // Node 2
                0.0, 0.0, 0.0;                                    // Node 3

            matrix x(10, 3);
            x.row(0) = x_vertices.row(0);
            x.row(1) = x_vertices.row(1);
            x.row(2) = x_vertices.row(2);
            x.row(3) = x_vertices.row(3);
            x.row(4) = 0.5 * (x_vertices.row(0) + x_vertices.row(1));
            x.row(5) = 0.5 * (x_vertices.row(1) + x_vertices.row(2));
            x.row(6) = 0.5 * (x_vertices.row(2) + x_vertices.row(3));
            x.row(7) = 0.5 * (x_vertices.row(3) + x_vertices.row(0));
            x.row(8) = 0.5 * (x_vertices.row(2) + x_vertices.row(0));
            x.row(9) = 0.5 * (x_vertices.row(3) + x_vertices.row(1));

            x.transposeInPlace();

            double const vol = tet10.quadrature().integrate(0.0,
                                                            [&](auto const& femval, auto const& l) {
                                                                auto const& [N, dN] = femval;

                                                                matrix3 const Jacobian = x * dN;

                                                                return Jacobian.determinant();
                                                            });

            REQUIRE(vol == Approx(1.0 / (6.0 * std::sqrt(2.0))));
        }
    }
    SECTION("Virtual methods check")
    {
        REQUIRE(make_volume_interpolation(element_topology::tetrahedron4, full())->nodes() == 4);
        REQUIRE(make_volume_interpolation(element_topology::tetrahedron10, full())->nodes() == 10);
    }
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
    SECTION("Six node - one point evaluation")
    {
        prism6 element(prism_quadrature::point::one);

        REQUIRE(element.nodes() == 6);
        REQUIRE(element.quadrature().points() == 1);

        element.quadrature().for_each([&](auto const& femval, auto const& l) {
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

        REQUIRE(element.nodes() == 6);
        REQUIRE(element.quadrature().points() == 6);

        element.quadrature().for_each([&](auto const& femval, auto const& l) {
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

        REQUIRE(element.nodes() == 15);
        REQUIRE(element.quadrature().points() == 6);

        element.quadrature().for_each([&](auto const& femval, auto const& l) {
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
        REQUIRE(make_volume_interpolation(element_topology::prism6, full())->nodes() == 6);
        REQUIRE(make_volume_interpolation(element_topology::prism15, full())->nodes() == 15);
    }
}
TEST_CASE("Unit sphere quadrature scheme test", "[unit_sphere_quadrature]")
{
    SECTION("BO21 rule test")
    {
        unit_sphere_quadrature unit_sphere(unit_sphere_quadrature::point::BO21);

        REQUIRE(unit_sphere.points() == 21);
        REQUIRE(ranges::accumulate(unit_sphere.weights(), 0.0) == Approx(1.0));

        for (auto const& coordinate : unit_sphere.coordinates())
        {
            auto const& [l, x, y, z] = coordinate;
            vector3 norm_check(x, y, z);
            REQUIRE(norm_check.norm() == Approx(1.0));
        }
    }
    SECTION("BO33 rule test")
    {
        unit_sphere_quadrature unit_sphere(unit_sphere_quadrature::point::BO33);

        REQUIRE(unit_sphere.points() == 33);
        REQUIRE(ranges::accumulate(unit_sphere.weights(), 0.0) == Approx(1.0));

        for (auto const& coordinate : unit_sphere.coordinates())
        {
            auto const& [l, x, y, z] = coordinate;
            vector3 norm_check(x, y, z);
            REQUIRE(norm_check.norm() == Approx(1.0));
        }
    }
    SECTION("FM900 rule test")
    {
        unit_sphere_quadrature unit_sphere(unit_sphere_quadrature::point::FM900);

        REQUIRE(unit_sphere.points() == 900);
        REQUIRE(ranges::accumulate(unit_sphere.weights(), 0.0) == Approx(12.5663706143));

        for (auto const& coordinate : unit_sphere.coordinates())
        {
            auto const& [l, x, y, z] = coordinate;
            vector3 norm_check(x, y, z);
            REQUIRE(norm_check.norm() == Approx(1.0));
        }
    }
}
