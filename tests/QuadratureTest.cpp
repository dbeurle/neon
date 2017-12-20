
#define CATCH_CONFIG_MAIN
#include <catch.hpp>

#include "quadrature/HexahedronQuadrature.hpp"
#include "quadrature/LineQuadrature.hpp"
#include "quadrature/PrismQuadrature.hpp"
#include "quadrature/QuadrilateralQuadrature.hpp"
#include "quadrature/TetrahedronQuadrature.hpp"
#include "quadrature/TriangleQuadrature.hpp"
#include "quadrature/UnitSphereQuadrature.hpp"

#include "interpolations/Hexahedron.hpp"
#include "interpolations/Line.hpp"
#include "interpolations/Quadrilateral.hpp"

#include "interpolations/Tetrahedron10.hpp"
#include "interpolations/Tetrahedron4.hpp"
#include "interpolations/Triangle3.hpp"
#include "interpolations/Triangle6.hpp"

#include <range/v3/numeric/accumulate.hpp>

using namespace neon;

constexpr auto ZERO_MARGIN = 1.0e-5;

TEST_CASE("Line quadrature scheme test", "[LineQuadrature]")
{
    SECTION("Line Gauss Quadrature")
    {
        // Check 1 and 8 point rule
        LineQuadrature l1(LineQuadrature::Rule::OnePoint);
        LineQuadrature l2(LineQuadrature::Rule::TwoPoint);
        LineQuadrature l3(LineQuadrature::Rule::ThreePoint);

        REQUIRE(l1.points() == 1);
        REQUIRE(l2.points() == 2);
        REQUIRE(l3.points() == 3);

        REQUIRE(ranges::accumulate(l1.weights(), 0.0) == Approx(2.0));
        REQUIRE(ranges::accumulate(l2.weights(), 0.0) == Approx(2.0));
        REQUIRE(ranges::accumulate(l3.weights(), 0.0) == Approx(2.0));
    }
}
TEST_CASE("Quadrilateral quadrature scheme test", "[QuadrilateralQuadrature]")
{
    SECTION("Quadrilateral Gauss Quadrature")
    {
        // Check 1 and 8 point rule
        QuadrilateralQuadrature q1(QuadrilateralQuadrature::Rule::OnePoint);
        QuadrilateralQuadrature q4(QuadrilateralQuadrature::Rule::FourPoint);
        QuadrilateralQuadrature q9(QuadrilateralQuadrature::Rule::NinePoint);

        REQUIRE(q1.points() == 1);
        REQUIRE(q4.points() == 4);
        REQUIRE(q9.points() == 9);

        REQUIRE(ranges::accumulate(q1.weights(), 0.0) == Approx(4.0));
        REQUIRE(ranges::accumulate(q4.weights(), 0.0) == Approx(4.0));
        REQUIRE(ranges::accumulate(q9.weights(), 0.0) == Approx(4.0));
    }
    SECTION("Quadrilateral4 interpolation function - one point")
    {
        Quadrilateral4 quad4(QuadrilateralQuadrature::Rule::OnePoint);

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
    SECTION("Quadrilateral4 surface area - one point")
    {
        Quadrilateral4 quad4(QuadrilateralQuadrature::Rule::OnePoint);

        matrix x(3, 4);
        x << 0.0, 1.0, 1.0, 0.0, //
            0.0, 0.0, 1.0, 1.0,  //
            0.0, 0.0, 0.0, 0.0;

        REQUIRE(quad4.compute_measure(x) == Approx(1.0));
    }
    SECTION("Quadrilateral4 interpolation function - four point")
    {
        Quadrilateral4 quad4(QuadrilateralQuadrature::Rule::FourPoint);

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
    SECTION("Quadrilateral4 surface area - four point")
    {
        Quadrilateral4 quad4(QuadrilateralQuadrature::Rule::FourPoint);

        matrix x(3, 4);
        x << 0.0, 1.0, 1.0, 0.0, //
            0.0, 0.0, 1.0, 1.0,  //
            0.0, 0.0, 0.0, 0.0;

        REQUIRE(quad4.compute_measure(x) == Approx(1.0));
    }
    SECTION("Quadrilateral8 interpolation function - four point")
    {
        Quadrilateral8 quad8(QuadrilateralQuadrature::Rule::FourPoint);

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
    SECTION("Quadrilateral8 surface area - four point")
    {
        Quadrilateral8 quad8(QuadrilateralQuadrature::Rule::FourPoint);

        matrix x(3, 8);
        x << 0.0, 1.0, 1.0, 0.0, 0.5, 1.0, 0.5, 0.0, //
            0.0, 0.0, 1.0, 1.0, 0.0, 0.5, 1.0, 0.5,  //
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0;

        REQUIRE(quad8.compute_measure(x) == Approx(1.0));
    }
    SECTION("Quadrilateral8 interpolation function - nine point")
    {
        Quadrilateral8 quad8(QuadrilateralQuadrature::Rule::NinePoint);

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
    SECTION("Quadrilateral8 interpolation function - nine point")
    {
        Quadrilateral8 quad8(QuadrilateralQuadrature::Rule::NinePoint);

        matrix x(3, 8);
        x << 0.0, 1.0, 1.0, 0.0, 0.5, 1.0, 0.5, 0.0, //
            0.0, 0.0, 1.0, 1.0, 0.0, 0.5, 1.0, 0.5,  //
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0;

        REQUIRE(quad8.compute_measure(x) == Approx(1.0));
    }
    SECTION("Quadrilateral9 interpolation function - four point")
    {
        Quadrilateral9 quad9(QuadrilateralQuadrature::Rule::FourPoint);

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
    SECTION("Quadrilateral9 surface area - four point")
    {
        Quadrilateral9 quad9(QuadrilateralQuadrature::Rule::FourPoint);

        matrix x(3, 9);
        x << 0.0, 1.0, 1.0, 0.0, 0.5, 1.0, 0.5, 0.0, 0.5, //
            0.0, 0.0, 1.0, 1.0, 0.0, 0.5, 1.0, 0.5, 0.5,  //
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0;

        REQUIRE(quad9.compute_measure(x) == Approx(1.0));
    }
    SECTION("Quadrilateral9 interpolation function - nine point")
    {
        Quadrilateral9 quad9(QuadrilateralQuadrature::Rule::NinePoint);

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
    SECTION("Quadrilateral9 interpolation function - nine point")
    {
        Quadrilateral9 quad9(QuadrilateralQuadrature::Rule::NinePoint);

        matrix x(3, 9);
        x << 0.0, 1.0, 1.0, 0.0, 0.5, 1.0, 0.5, 0.0, 0.5, //
            0.0, 0.0, 1.0, 1.0, 0.0, 0.5, 1.0, 0.5, 0.5,  //
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0;

        REQUIRE(quad9.compute_measure(x) == Approx(1.0));
    }
}
TEST_CASE("Triangle quadrature scheme test", "[TriangleQuadrature]")
{
    SECTION("Triangle Gauss Quadrature")
    {
        // Check 1 and 8 point rule
        TriangleQuadrature t1(TriangleQuadrature::Rule::OnePoint);
        TriangleQuadrature t3(TriangleQuadrature::Rule::ThreePoint);
        TriangleQuadrature t4(TriangleQuadrature::Rule::FourPoint);

        REQUIRE(t1.points() == 1);
        REQUIRE(t3.points() == 3);
        REQUIRE(t4.points() == 4);

        REQUIRE(ranges::accumulate(t1.weights(), 0.0) == Approx(0.5));
        REQUIRE(ranges::accumulate(t3.weights(), 0.0) == Approx(0.5));
        REQUIRE(ranges::accumulate(t4.weights(), 0.0) == Approx(0.5));
    }
    SECTION("Triangle3 interpolation function - one point")
    {
        Triangle3 tri3(TriangleQuadrature::Rule::OnePoint);

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
    SECTION("Triangle6 interpolation function - three point")
    {
        Triangle6 tri6(TriangleQuadrature::Rule::ThreePoint);

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
    SECTION("Triangle6 interpolation function - four point")
    {
        Triangle6 tri6(TriangleQuadrature::Rule::FourPoint);

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
    SECTION("Triangle3 surface area")
    {
        Triangle3 tri3(TriangleQuadrature::Rule::OnePoint);

        matrix x(3, 3);
        x << 0.0, 0.0, 0.0, //
            1.0, 0.0, 0.0,  //
            0.0, 1.0, 0.0;
        x.transposeInPlace();

        REQUIRE(tri3.compute_measure(x) == Approx(0.5));
    }
}
TEST_CASE("Hexahedron quadrature scheme test", "[HexahedronQuadrature]")
{
    SECTION("Hexahedron Gauss Quadrature")
    {
        // Check 1 and 8 point rule
        HexahedronQuadrature One(HexahedronQuadrature::Rule::OnePoint);
        HexahedronQuadrature Six(HexahedronQuadrature::Rule::SixPoint);
        HexahedronQuadrature Eight(HexahedronQuadrature::Rule::EightPoint);
        HexahedronQuadrature TwentySeven(HexahedronQuadrature::Rule::TwentySevenPoint);

        REQUIRE(One.points() == 1);
        REQUIRE(Six.points() == 6);
        REQUIRE(Eight.points() == 8);
        REQUIRE(TwentySeven.points() == 27);

        REQUIRE(ranges::accumulate(One.weights(), 0.0) == Approx(8.0));
        REQUIRE(ranges::accumulate(Six.weights(), 0.0) == Approx(8.0));
        REQUIRE(ranges::accumulate(Eight.weights(), 0.0) == Approx(8.0));
        REQUIRE(ranges::accumulate(TwentySeven.weights(), 0.0) == Approx(8.0));
    }
    SECTION("Hexahedron8 OnePoint Evaluation")
    {
        Hexahedron8 hex8(HexahedronQuadrature::Rule::OnePoint);

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
    SECTION("Hexahedron8 EightPoint Evaluation")
    {
        Hexahedron8 hex8(HexahedronQuadrature::Rule::EightPoint);

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
    SECTION("Hexahedron8 volume evaluation")
    {
        Hexahedron8 hex8(HexahedronQuadrature::Rule::EightPoint);

        matrix x(3, 8);
        x << 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, //
            0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0,  //
            0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0;

        REQUIRE(hex8.compute_measure(x) == Approx(1.0));
    }
    SECTION("Hexahedron20 OnePoint Evaluation")
    {
        Hexahedron20 hex20(HexahedronQuadrature::Rule::OnePoint);

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
    SECTION("Hexahedron20 SixPoint Evaluation")
    {
        Hexahedron20 hex20(HexahedronQuadrature::Rule::SixPoint);

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
    SECTION("Hexahedron20 volume evaluation")
    {
        Hexahedron20 hex20(HexahedronQuadrature::Rule::EightPoint);

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
    SECTION("Hexahedron20 EightPoint Evaluation")
    {
        Hexahedron20 hex20(HexahedronQuadrature::Rule::EightPoint);

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
    SECTION("Hexahedron20 TwentySevenPoint Evaluation")
    {
        Hexahedron20 hex20(HexahedronQuadrature::Rule::TwentySevenPoint);

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
    SECTION("Hexahedron27 OnePoint Evaluation")
    {
        Hexahedron27 hex27(HexahedronQuadrature::Rule::OnePoint);

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
    SECTION("Hexahedron27 SixPoint Evaluation")
    {
        Hexahedron27 hex27(HexahedronQuadrature::Rule::SixPoint);

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
    SECTION("Hexahedron27 EightPoint Evaluation")
    {
        Hexahedron27 hex27(HexahedronQuadrature::Rule::EightPoint);

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
    SECTION("Hexahedron27 TwentySevenPoint Evaluation")
    {
        Hexahedron27 hex27(HexahedronQuadrature::Rule::TwentySevenPoint);

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
    SECTION("Hexahedron27 volume evaluation")
    {
        SECTION("Six point rule")
        {
            Hexahedron27 hex27(HexahedronQuadrature::Rule::SixPoint);

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
            Hexahedron27 hex27(HexahedronQuadrature::Rule::EightPoint);

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
TEST_CASE("Tetrahedron quadrature scheme test", "[TetrahedronQuadrature]")
{
    SECTION("Tetrahedron Gauss Quadrature")
    {
        TetrahedronQuadrature One(TetrahedronQuadrature::Rule::OnePoint);
        TetrahedronQuadrature Four(TetrahedronQuadrature::Rule::FourPoint);
        TetrahedronQuadrature Five(TetrahedronQuadrature::Rule::FivePoint);

        // Check the number of quadrature points are correctly set
        REQUIRE(One.points() == 1);
        REQUIRE(Four.points() == 4);
        REQUIRE(Five.points() == 5);

        // Check the weightings add to 1/6 for each scheme
        REQUIRE(ranges::accumulate(One.weights(), 0.0) == Approx(1.0 / 6.0));
        REQUIRE(ranges::accumulate(Four.weights(), 0.0) == Approx(1.0 / 6.0));
        REQUIRE(ranges::accumulate(Five.weights(), 0.0) == Approx(1.0 / 6.0));
    }
    SECTION("Tetrahedron4 OnePoint Evaluation")
    {
        Tetrahedron4 tet4(TetrahedronQuadrature::Rule::OnePoint);

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
    SECTION("Tetrahedron10 OnePoint Evaluation")
    {
        Tetrahedron10 tet10(TetrahedronQuadrature::Rule::OnePoint);

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
    SECTION("Tetrahedron10 FourPoint Evaluation")
    {
        Tetrahedron10 tet10(TetrahedronQuadrature::Rule::FourPoint);

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
    SECTION("Tetrahedron10 FivePoint Evaluation")
    {
        Tetrahedron10 tet10(TetrahedronQuadrature::Rule::FivePoint);

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
    SECTION("Tetrahedron10 Jacobian and volume check")
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
                Tetrahedron10 tet10(TetrahedronQuadrature::Rule::OnePoint);

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
                Tetrahedron10 tet10(TetrahedronQuadrature::Rule::FourPoint);

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
                Tetrahedron10 tet10(TetrahedronQuadrature::Rule::FivePoint);

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
            Tetrahedron10 tet10(TetrahedronQuadrature::Rule::FourPoint);

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
}
TEST_CASE("Prism quadrature scheme test", "[PrismQuadrature]")
{
    SECTION("Prism Gauss Quadrature")
    {
        // Check 1 and 6 point rule
        PrismQuadrature p1(PrismQuadrature::Rule::OnePoint);
        PrismQuadrature p6(PrismQuadrature::Rule::SixPoint);

        REQUIRE(p1.points() == 1);
        REQUIRE(p6.points() == 6);

        REQUIRE(ranges::accumulate(p1.weights(), 0.0) == Approx(4.0));
        REQUIRE(ranges::accumulate(p6.weights(), 0.0) == Approx(6.0));
    }
}
TEST_CASE("Unit sphere quadrature scheme test", "[UnitSphereQuadrature]")
{
    SECTION("BO21 rule test")
    {
        UnitSphereQuadrature unit_sphere(UnitSphereQuadrature::Rule::BO21);

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
        UnitSphereQuadrature unit_sphere(UnitSphereQuadrature::Rule::BO33);

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
        UnitSphereQuadrature unit_sphere(UnitSphereQuadrature::Rule::FM900);

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
