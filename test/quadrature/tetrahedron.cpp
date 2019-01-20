
#include <catch2/catch.hpp>

#include "quadrature/tetrahedron/jinyun.hpp"
#include "quadrature/tetrahedron/witherden_vincent.hpp"

#include "interpolations/tetrahedron.hpp"

#include "interpolations/interpolation_factory.hpp"
#include "mesh/element_topology.hpp"

#include "io/json.hpp"

#include <numeric>

using namespace neon;

json full() { return json::parse("{\"element_options\" : {\"quadrature\" : \"full\"}}"); }
json reduced() { return json::parse("{\"element_options\" : {\"quadrature\" : \"reduced\"}}"); }

constexpr auto ZERO_MARGIN = 1.0e-5;

TEST_CASE("tetrahedron quadrature")
{
    using namespace neon::quadrature::tetrahedron;

    SECTION("jinyun")
    {
        // tetrahedron_quadrature One(1);
        // tetrahedron_quadrature Four(2);
        // tetrahedron_quadrature Five(3);
        //
        // // Check the number of quadrature points are correctly set
        // REQUIRE(One.points() == 1);
        // REQUIRE(Four.points() == 4);
        // REQUIRE(Five.points() == 5);
        //
        // // Check the weightings add to 1/6 for each scheme
        // REQUIRE(ranges::accumulate(One.weights(), 0.0) == Approx(1.0 / 6.0));
        // REQUIRE(ranges::accumulate(Four.weights(), 0.0) == Approx(1.0 / 6.0));
        // REQUIRE(ranges::accumulate(Five.weights(), 0.0) == Approx(1.0 / 6.0));
    }
    SECTION("witherden_vincent")
    {
        // tetrahedron_quadrature One(1);
        // tetrahedron_quadrature Four(2);
        // tetrahedron_quadrature Five(3);
        //
        // // Check the number of quadrature points are correctly set
        // REQUIRE(One.points() == 1);
        // REQUIRE(Four.points() == 4);
        // REQUIRE(Five.points() == 5);
        //
        // // Check the weightings add to 1/6 for each scheme
        // REQUIRE(ranges::accumulate(One.weights(), 0.0) == Approx(1.0 / 6.0));
        // REQUIRE(ranges::accumulate(Four.weights(), 0.0) == Approx(1.0 / 6.0));
        // REQUIRE(ranges::accumulate(Five.weights(), 0.0) == Approx(1.0 / 6.0));
    }
}
TEST_CASE("tetrahedron shape functions")
{
    using namespace neon::quadrature::tetrahedron;

    SECTION("tetrahedron4 jinyun quadrature")
    {
        tetrahedron4 element;

        jinyun scheme(1);

        REQUIRE(element.number_of_nodes() == 4);

        for (auto const& coordinate : scheme.coordinates())
        {
            auto const [N, dN] = element.evaluate(coordinate);

            REQUIRE(N.sum() == Approx(1.0));

            REQUIRE(dN.col(0).sum() == Approx(0.0).margin(ZERO_MARGIN));
            REQUIRE(dN.col(1).sum() == Approx(0.0).margin(ZERO_MARGIN));
            REQUIRE(dN.col(2).sum() == Approx(0.0).margin(ZERO_MARGIN));
        }
    }
    SECTION("tetrahedron4 witherden vincent quadrature")
    {
        tetrahedron4 element;

        witherden_vincent scheme(1);

        REQUIRE(element.number_of_nodes() == 4);

        for (auto const& coordinate : scheme.coordinates())
        {
            auto const [N, dN] = element.evaluate(coordinate);

            REQUIRE(N.sum() == Approx(1.0));

            REQUIRE(dN.col(0).sum() == Approx(0.0).margin(ZERO_MARGIN));
            REQUIRE(dN.col(1).sum() == Approx(0.0).margin(ZERO_MARGIN));
            REQUIRE(dN.col(2).sum() == Approx(0.0).margin(ZERO_MARGIN));
        }
    }
    SECTION("tetrahedron10 jinyun quadrature")
    {
        tetrahedron10 element;

        jinyun scheme(2);

        REQUIRE(element.number_of_nodes() == 10);

        for (auto const& coordinate : scheme.coordinates())
        {
            auto const [N, dN] = element.evaluate(coordinate);

            REQUIRE(N.sum() == Approx(1.0));

            REQUIRE(dN.col(0).sum() == Approx(0.0).margin(ZERO_MARGIN));
            REQUIRE(dN.col(1).sum() == Approx(0.0).margin(ZERO_MARGIN));
            REQUIRE(dN.col(2).sum() == Approx(0.0).margin(ZERO_MARGIN));
        }
    }
    SECTION("tetrahedron10 witherden vincent quadrature")
    {
        tetrahedron10 element;

        witherden_vincent scheme(2);

        REQUIRE(element.number_of_nodes() == 10);

        for (auto const& coordinate : scheme.coordinates())
        {
            auto const [N, dN] = element.evaluate(coordinate);

            REQUIRE(N.sum() == Approx(1.0));

            REQUIRE(dN.col(0).sum() == Approx(0.0).margin(ZERO_MARGIN));
            REQUIRE(dN.col(1).sum() == Approx(0.0).margin(ZERO_MARGIN));
            REQUIRE(dN.col(2).sum() == Approx(0.0).margin(ZERO_MARGIN));
        }
    }

    // SECTION("tetrahedron10 Jacobian and volume check")
    // {
    //     SECTION("Parent element")
    //     {
    //         // Check that we have a Jacobian of one with a parent element
    //         matrix x(10, 3);
    //         x << 1.0, 0.0, 0.0, //
    //             0.0, 1.0, 0.0,  //
    //             0.0, 0.0, 1.0,  //
    //             0.0, 0.0, 0.0,  //
    //             0.5, 0.5, 0.0,  //
    //             0.0, 0.5, 0.5,  //
    //             0.0, 0.0, 0.5,  //
    //             0.5, 0.0, 0.0,  //
    //             0.5, 0.0, 0.5,  //
    //             0.0, 0.5, 0.0;
    //         x.transposeInPlace();
    //
    //         SECTION("One point")
    //         {
    //             tetrahedron10 tet10(tetrahedron_quadrature::point::one);
    //
    //             tet10.quadrature().for_each([&](auto const& value, auto) {
    //                 auto const [N, dN] = value;
    //
    //                 matrix3 const Jacobian = x * dN;
    //
    //                 REQUIRE(Jacobian.determinant() == Approx(1.0));
    //             });
    //             double vol = tet10.quadrature().integrate(0.0, [&](auto const& value, auto) {
    //                 auto const [N, dN] = value;
    //
    //                 matrix3 const Jacobian = x * dN;
    //
    //                 return Jacobian.determinant();
    //             });
    //             REQUIRE(vol == Approx(1.0 / 6.0));
    //         }
    //         SECTION("Four point")
    //         {
    //             tetrahedron10 tet10(tetrahedron_quadrature::point::four);
    //
    //             tet10.quadrature().for_each([&](auto const& value, auto) {
    //                 auto const [N, dN] = value;
    //
    //                 matrix3 const Jacobian = x * dN;
    //
    //                 REQUIRE(Jacobian.determinant() == Approx(1.0));
    //             });
    //
    //             double vol = tet10.quadrature().integrate(0.0, [&](auto const& value, auto) {
    //                 auto const [N, dN] = value;
    //
    //                 matrix3 const Jacobian = x * dN;
    //
    //                 return Jacobian.determinant();
    //             });
    //             REQUIRE(vol == Approx(1.0 / 6.0));
    //         }
    //         SECTION("Five point")
    //         {
    //             tetrahedron10 tet10(tetrahedron_quadrature::point::five);
    //
    //             tet10.quadrature().for_each([&](auto const& value, auto) {
    //                 auto const [N, dN] = value;
    //
    //                 matrix3 const Jacobian = x * dN;
    //
    //                 REQUIRE(Jacobian.determinant() == Approx(1.0));
    //             });
    //
    //             double vol = tet10.quadrature().integrate(0.0, [&](auto const& value, auto) {
    //                 auto const [N, dN] = value;
    //
    //                 matrix3 const Jacobian = x * dN;
    //
    //                 return Jacobian.determinant();
    //             });
    //             REQUIRE(vol == Approx(1.0 / 6.0));
    //         }
    //     }
    //     SECTION("Volume of tetrahedron with unit length sides")
    //     {
    //         tetrahedron10 tet10(tetrahedron_quadrature::point::four);
    //
    //         matrix x_vertices(4, 3);
    //         x_vertices << 1.0, 0.0, 0.0,                          // Node 0
    //             0.5, 0.5 * std::sqrt(3.0), 0.0,                   // Node 1
    //             0.5, std::sqrt(1.0 / 12.0), std::sqrt(2.0 / 3.0), // Node 2
    //             0.0, 0.0, 0.0;                                    // Node 3
    //
    //         matrix x(10, 3);
    //         x.row(0) = x_vertices.row(0);
    //         x.row(1) = x_vertices.row(1);
    //         x.row(2) = x_vertices.row(2);
    //         x.row(3) = x_vertices.row(3);
    //         x.row(4) = 0.5 * (x_vertices.row(0) + x_vertices.row(1));
    //         x.row(5) = 0.5 * (x_vertices.row(1) + x_vertices.row(2));
    //         x.row(6) = 0.5 * (x_vertices.row(2) + x_vertices.row(3));
    //         x.row(7) = 0.5 * (x_vertices.row(3) + x_vertices.row(0));
    //         x.row(8) = 0.5 * (x_vertices.row(2) + x_vertices.row(0));
    //         x.row(9) = 0.5 * (x_vertices.row(3) + x_vertices.row(1));
    //
    //         x.transposeInPlace();
    //
    //         double const vol = tet10.quadrature().integrate(0.0, [&](auto const& value, auto) {
    //             auto const [N, dN] = value;
    //
    //             matrix3 const Jacobian = x * dN;
    //
    //             return Jacobian.determinant();
    //         });
    //
    //         REQUIRE(vol == Approx(1.0 / (6.0 * std::sqrt(2.0))));
    //     }
    // }
    // SECTION("Virtual methods check")
    // {
    //     REQUIRE(make_volume_interpolation(element_topology::tetrahedron4, full())->number_of_nodes()
    //             == 4);
    //     REQUIRE(make_volume_interpolation(element_topology::tetrahedron10, full())->number_of_nodes()
    //             == 10);
    // }
}

// REQUIRE(tet10.local_quadrature_extrapolation().rows() == 10);
// REQUIRE(tet10.local_quadrature_extrapolation().cols() == 4);
// REQUIRE(tet10.local_quadrature_extrapolation().allFinite());
