
#include <catch2/catch.hpp>

#include "quadrature/tetrahedron/tetrahedron_quadrature.hpp"
#include "interpolations/tetrahedron.hpp"

#include "interpolations/interpolation_factory.hpp"
#include "mesh/element_topology.hpp"

#include "io/json.hpp"

#include <range/v3/numeric/accumulate.hpp>

using namespace neon;

json full() { return json::parse("{\"element_options\" : {\"quadrature\" : \"full\"}}"); }
json reduced() { return json::parse("{\"element_options\" : {\"quadrature\" : \"reduced\"}}"); }

constexpr auto ZERO_MARGIN = 1.0e-5;

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

        REQUIRE(tet4.number_of_nodes() == 4);
        REQUIRE(tet4.quadrature().points() == 1);

        tet4.quadrature().for_each([&](auto const& femval, auto) {
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

        REQUIRE(tet10.number_of_nodes() == 10);

        REQUIRE(tet10.quadrature().points() == 1);

        tet10.quadrature().for_each([&](auto const& femval, auto) {
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

        REQUIRE(tet10.number_of_nodes() == 10);

        REQUIRE(tet10.quadrature().points() == 4);

        tet10.quadrature().for_each([&](auto const& femval, auto) {
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

        REQUIRE(tet10.number_of_nodes() == 10);

        REQUIRE(tet10.quadrature().points() == 5);

        tet10.quadrature().for_each([&](auto const& femval, auto) {
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

                tet10.quadrature().for_each([&](auto const& femval, auto) {
                    auto const& [N, rhea] = femval;

                    matrix3 const Jacobian = x * rhea;

                    REQUIRE(Jacobian.determinant() == Approx(1.0));
                });
                double vol = tet10.quadrature().integrate(0.0, [&](auto const& femval, auto) {
                    auto const& [N, rhea] = femval;

                    matrix3 const Jacobian = x * rhea;

                    return Jacobian.determinant();
                });
                REQUIRE(vol == Approx(1.0 / 6.0));
            }
            SECTION("Four point")
            {
                tetrahedron10 tet10(tetrahedron_quadrature::point::four);

                tet10.quadrature().for_each([&](auto const& femval, auto) {
                    auto const& [N, rhea] = femval;

                    matrix3 const Jacobian = x * rhea;

                    REQUIRE(Jacobian.determinant() == Approx(1.0));
                });

                double vol = tet10.quadrature().integrate(0.0, [&](auto const& femval, auto) {
                    auto const& [N, rhea] = femval;

                    matrix3 const Jacobian = x * rhea;

                    return Jacobian.determinant();
                });
                REQUIRE(vol == Approx(1.0 / 6.0));
            }
            SECTION("Five point")
            {
                tetrahedron10 tet10(tetrahedron_quadrature::point::five);

                tet10.quadrature().for_each([&](auto const& femval, auto) {
                    auto const& [N, rhea] = femval;

                    matrix3 const Jacobian = x * rhea;

                    REQUIRE(Jacobian.determinant() == Approx(1.0));
                });

                double vol = tet10.quadrature().integrate(0.0, [&](auto const& femval, auto) {
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

            double const vol = tet10.quadrature().integrate(0.0, [&](auto const& femval, auto) {
                auto const& [N, dN] = femval;

                matrix3 const Jacobian = x * dN;

                return Jacobian.determinant();
            });

            REQUIRE(vol == Approx(1.0 / (6.0 * std::sqrt(2.0))));
        }
    }
    SECTION("Virtual methods check")
    {
        REQUIRE(make_volume_interpolation(element_topology::tetrahedron4, full())->number_of_nodes()
                == 4);
        REQUIRE(make_volume_interpolation(element_topology::tetrahedron10, full())->number_of_nodes()
                == 10);
    }
}
