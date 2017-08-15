
#define CATCH_CONFIG_MAIN // This tells Catch to provide a main() - only do this in one
                          // cpp file
#include <catch.hpp>

#include "quadrature/HexahedronQuadrature.hpp"
#include "quadrature/PrismQuadrature.hpp"
#include "quadrature/QuadrilateralQuadrature.hpp"
#include "quadrature/TetrahedronQuadrature.hpp"
#include "quadrature/TriangleQuadrature.hpp"
#include "quadrature/UnitSphereQuadrature.hpp"

#include "interpolations/Hexahedron8.hpp"
#include "interpolations/Quadrilateral4.hpp"
#include "interpolations/Quadrilateral8.hpp"
#include "interpolations/Tetrahedron10.hpp"
#include "interpolations/Tetrahedron4.hpp"
#include "interpolations/Triangle3.hpp"
#include "interpolations/Triangle6.hpp"

#include <range/v3/numeric.hpp>

using namespace neon;

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

            auto const & [ N, rhea ] = femval;

            REQUIRE(N.sum() == Approx(1.0));

            REQUIRE(rhea.col(0).sum() == Approx(0.0));
            REQUIRE(rhea.col(1).sum() == Approx(0.0));
        });
    }
    SECTION("Quadrilateral4 interpolation function - four point")
    {
        Quadrilateral4 quad4(QuadrilateralQuadrature::Rule::FourPoint);

        REQUIRE(quad4.nodes() == 4);
        REQUIRE(quad4.quadrature().points() == 4);
        REQUIRE(quad4.local_quadrature_extrapolation().rows() == 4);
        REQUIRE(quad4.local_quadrature_extrapolation().cols() == 4);

        quad4.quadrature().for_each([&](auto const& femval, auto const& l) {

            auto const & [ N, rhea ] = femval;

            REQUIRE(N.sum() == Approx(1.0));

            REQUIRE(rhea.col(0).sum() == Approx(0.0));
            REQUIRE(rhea.col(1).sum() == Approx(0.0));
        });
    }
    SECTION("Quadrilateral8 interpolation function - four point")
    {
        Quadrilateral8 quad8(QuadrilateralQuadrature::Rule::FourPoint);

        REQUIRE(quad8.nodes() == 8);
        REQUIRE(quad8.quadrature().points() == 4);
        REQUIRE(quad8.local_quadrature_extrapolation().rows() == 8);
        REQUIRE(quad8.local_quadrature_extrapolation().cols() == 4);

        quad8.quadrature().for_each([&](auto const& femval, auto const& l) {

            auto const & [ N, rhea ] = femval;

            REQUIRE(N.sum() == Approx(1.0));

            REQUIRE(rhea.col(0).sum() == Approx(0.0));
            REQUIRE(rhea.col(1).sum() == Approx(0.0));
        });
    }
    SECTION("Quadrilateral8 interpolation function - nine point")
    {
        Quadrilateral8 quad8(QuadrilateralQuadrature::Rule::NinePoint);

        REQUIRE(quad8.nodes() == 8);
        REQUIRE(quad8.quadrature().points() == 9);
        REQUIRE(quad8.local_quadrature_extrapolation().rows() == 8);
        REQUIRE(quad8.local_quadrature_extrapolation().cols() == 9);

        quad8.quadrature().for_each([&](auto const& femval, auto const& l) {

            auto const & [ N, rhea ] = femval;

            REQUIRE(N.sum() == Approx(1.0));

            REQUIRE(rhea.col(0).sum() == Approx(0.0));
            REQUIRE(rhea.col(1).sum() == Approx(0.0));
        });
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

            auto const & [ N, rhea ] = femval;

            REQUIRE(N.sum() == Approx(1.0));

            REQUIRE(rhea.col(0).sum() == Approx(0.0));
            REQUIRE(rhea.col(1).sum() == Approx(0.0));
        });
    }
    SECTION("Triangle6 interpolation function - three point")
    {
        Triangle6 tri6(TriangleQuadrature::Rule::ThreePoint);

        REQUIRE(tri6.nodes() == 6);
        REQUIRE(tri6.quadrature().points() == 3);

        tri6.quadrature().for_each([&](auto const& femval, auto const& l) {

            auto const & [ N, rhea ] = femval;

            REQUIRE(N.sum() == Approx(1.0));

            REQUIRE(rhea.col(0).sum() == Approx(0.0));
            REQUIRE(rhea.col(1).sum() == Approx(0.0));
        });
    }
    SECTION("Triangle6 interpolation function - four point")
    {
        Triangle6 tri6(TriangleQuadrature::Rule::FourPoint);

        REQUIRE(tri6.nodes() == 6);
        REQUIRE(tri6.quadrature().points() == 4);

        tri6.quadrature().for_each([&](auto const& femval, auto const& l) {

            auto const & [ N, rhea ] = femval;

            REQUIRE(N.sum() == Approx(1.0));

            REQUIRE(rhea.col(0).sum() == Approx(0.0));
            REQUIRE(rhea.col(1).sum() == Approx(0.0));
        });
    }
}
TEST_CASE("Hexahedron quadrature scheme test", "[HexahedronQuadrature]")
{
    SECTION("Hexahedron Gauss Quadrature")
    {
        // Check 1 and 8 point rule
        HexahedronQuadrature One(HexahedronQuadrature::Rule::OnePoint);
        HexahedronQuadrature Eight(HexahedronQuadrature::Rule::EightPoint);

        REQUIRE(One.points() == 1);
        REQUIRE(Eight.points() == 8);

        REQUIRE(ranges::accumulate(One.weights(), 0.0) == Approx(8.0));
        REQUIRE(ranges::accumulate(Eight.weights(), 0.0) == Approx(8.0));
    }

    SECTION("Hexahedron8 OnePoint Evaluation")
    {
        Hexahedron8 hex8(HexahedronQuadrature::Rule::OnePoint);

        REQUIRE(hex8.nodes() == 8);

        // Check the quadrature is assigned correctly
        REQUIRE(hex8.quadrature().points() == 1);

        hex8.quadrature().for_each([&](auto const& femval, auto const& l) {

            auto const & [ N, rhea ] = femval;

            // Check every entry
            for (int i = 0; i < 8; i++) REQUIRE(N(i) == Approx(0.125));

            REQUIRE(N.sum() == Approx(1.0));

            REQUIRE(rhea.col(0).sum() == Approx(0.0));
            REQUIRE(rhea.col(1).sum() == Approx(0.0));
            REQUIRE(rhea.col(2).sum() == Approx(0.0));

        });

        REQUIRE(hex8.local_quadrature_extrapolation().rows() == 8);
        REQUIRE(hex8.local_quadrature_extrapolation().cols() == 1);
    }
    SECTION("Hexahedron8 EightPoint Evaluation")
    {
        Hexahedron8 hex8(HexahedronQuadrature::Rule::EightPoint);

        // Check the nodes are correct
        REQUIRE(hex8.nodes() == 8);

        // Check the quaduratre is assigned correctly
        REQUIRE(hex8.quadrature().points() == 8);

        hex8.quadrature().for_each([&](auto const& femval, auto const& l) {

            auto const & [ N, rhea ] = femval;

            REQUIRE(N.sum() == Approx(1.0));

            REQUIRE(rhea.col(0).sum() == Approx(0.0));
            REQUIRE(rhea.col(1).sum() == Approx(0.0));
            REQUIRE(rhea.col(2).sum() == Approx(0.0));

        });

        REQUIRE(hex8.local_quadrature_extrapolation().rows() == 8);
        REQUIRE(hex8.local_quadrature_extrapolation().cols() == 8);
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

            auto const & [ N, rhea ] = femval;

            REQUIRE(N.sum() == Approx(1.0));

            REQUIRE(rhea.col(0).sum() == Approx(0.0));
            REQUIRE(rhea.col(1).sum() == Approx(0.0));
            REQUIRE(rhea.col(2).sum() == Approx(0.0));
        });

        REQUIRE(tet4.local_quadrature_extrapolation().rows() == 4);
        REQUIRE(tet4.local_quadrature_extrapolation().cols() == 1);
    }
    SECTION("Tetrahedron10 OnePoint Evaluation")
    {
        Tetrahedron10 tet10(TetrahedronQuadrature::Rule::OnePoint);

        REQUIRE(tet10.nodes() == 10);

        REQUIRE(tet10.quadrature().points() == 1);

        tet10.quadrature().for_each([&](auto const& femval, auto const& l) {

            auto const & [ N, rhea ] = femval;

            REQUIRE(N.sum() == Approx(1.0));

            REQUIRE(rhea.col(0).sum() == Approx(0.0));
            REQUIRE(rhea.col(1).sum() == Approx(0.0));
            REQUIRE(rhea.col(2).sum() == Approx(0.0));
        });

        REQUIRE(tet10.local_quadrature_extrapolation().rows() == 10);
        REQUIRE(tet10.local_quadrature_extrapolation().cols() == 1);
    }
    SECTION("Tetrahedron10 FourPoint Evaluation")
    {
        Tetrahedron10 tet10(TetrahedronQuadrature::Rule::FourPoint);

        REQUIRE(tet10.nodes() == 10);

        REQUIRE(tet10.quadrature().points() == 4);

        tet10.quadrature().for_each([&](auto const& femval, auto const& l) {

            auto const & [ N, rhea ] = femval;

            REQUIRE(N.sum() == Approx(1.0));

            REQUIRE(rhea.col(0).sum() == Approx(0.0));
            REQUIRE(rhea.col(1).sum() == Approx(0.0));
            REQUIRE(rhea.col(2).sum() == Approx(0.0));
        });

        REQUIRE(tet10.local_quadrature_extrapolation().rows() == 10);
        REQUIRE(tet10.local_quadrature_extrapolation().cols() == 4);
    }
    SECTION("Tetrahedron10 FivePoint Evaluation")
    {
        Tetrahedron10 tet10(TetrahedronQuadrature::Rule::FivePoint);

        REQUIRE(tet10.nodes() == 10);

        REQUIRE(tet10.quadrature().points() == 5);

        tet10.quadrature().for_each([&](auto const& femval, auto const& l) {

            auto const & [ N, rhea ] = femval;

            REQUIRE(N.sum() == Approx(1.0));

            REQUIRE(rhea.col(0).sum() == Approx(0.0));
            REQUIRE(rhea.col(1).sum() == Approx(0.0));
            REQUIRE(rhea.col(2).sum() == Approx(0.0));
        });

        REQUIRE(tet10.local_quadrature_extrapolation().rows() == 10);
        REQUIRE(tet10.local_quadrature_extrapolation().cols() == 5);
    }
    SECTION("Tetrahedron10 Jacobian and volume check")
    {
        SECTION("Parent element")
        {
            // Check that we have a Jacobian of one with a parent element
            Matrix x(10, 3);
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

                    auto const & [ N, rhea ] = femval;

                    Matrix3 const Jacobian = x * rhea;

                    REQUIRE(Jacobian.determinant() == Approx(1.0));

                });
                double vol = tet10.quadrature()
                                 .integrate(0.0, [&](auto const& femval, auto const& l) {
                                     auto const & [ N, rhea ] = femval;

                                     Matrix3 const Jacobian = x * rhea;

                                     return Jacobian.determinant();
                                 });
                REQUIRE(vol == Approx(1.0 / 6.0));
            }
            SECTION("Four point")
            {
                Tetrahedron10 tet10(TetrahedronQuadrature::Rule::FourPoint);

                tet10.quadrature().for_each([&](auto const& femval, auto const& l) {

                    auto const & [ N, rhea ] = femval;

                    Matrix3 const Jacobian = x * rhea;

                    REQUIRE(Jacobian.determinant() == Approx(1.0));

                });

                double vol = tet10.quadrature()
                                 .integrate(0.0, [&](auto const& femval, auto const& l) {
                                     auto const & [ N, rhea ] = femval;

                                     Matrix3 const Jacobian = x * rhea;

                                     return Jacobian.determinant();
                                 });
                REQUIRE(vol == Approx(1.0 / 6.0));
            }
            SECTION("Five point")
            {
                Tetrahedron10 tet10(TetrahedronQuadrature::Rule::FivePoint);

                tet10.quadrature().for_each([&](auto const& femval, auto const& l) {

                    auto const & [ N, rhea ] = femval;

                    Matrix3 const Jacobian = x * rhea;

                    REQUIRE(Jacobian.determinant() == Approx(1.0));

                });

                double vol = tet10.quadrature()
                                 .integrate(0.0, [&](auto const& femval, auto const& l) {
                                     auto const & [ N, rhea ] = femval;

                                     Matrix3 const Jacobian = x * rhea;

                                     return Jacobian.determinant();
                                 });
                REQUIRE(vol == Approx(1.0 / 6.0));
            }
        }
        SECTION("Volume of tetrahedron with unit length sides")
        {
            Tetrahedron10 tet10(TetrahedronQuadrature::Rule::FourPoint);

            Matrix x_vertices(4, 3);
            x_vertices << 1.0, 0.0, 0.0,                          // Node 0
                0.5, 0.5 * std::sqrt(3.0), 0.0,                   // Node 1
                0.5, std::sqrt(1.0 / 12.0), std::sqrt(2.0 / 3.0), // Node 2
                0.0, 0.0, 0.0;                                    // Node 3

            Matrix x(10, 3);
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

            double const vol = tet10.quadrature()
                                   .integrate(0.0, [&](auto const& femval, auto const& l) {
                                       auto const & [ N, dN ] = femval;

                                       Matrix3 const Jacobian = x * dN;

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
    SECTION("Unit sphere unit test")
    {
        UnitSphereQuadrature unit_sphere;

        REQUIRE(unit_sphere.points() == 21);
        REQUIRE(ranges::accumulate(unit_sphere.weights(), 0.0) == Approx(1.0));
    }
}
