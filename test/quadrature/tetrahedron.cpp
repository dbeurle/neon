
#include <catch2/catch.hpp>

#include "quadrature/tetrahedron/jinyun.hpp"
#include "quadrature/tetrahedron/witherden_vincent.hpp"

#include "interpolations/tetrahedron.hpp"

#include "interpolations/interpolation_factory.hpp"
#include "mesh/element_topology.hpp"

#include <numeric>

using namespace neon;

constexpr auto ZERO_MARGIN = 1.0e-5;

template <typename ShapeFunction, typename Quadrature>
void check_shape_functions(ShapeFunction&& shape_function, Quadrature&& integration)
{
    for (auto const& coordinate : integration.coordinates())
    {
        auto const [N, dN] = shape_function.evaluate(coordinate);

        REQUIRE(N.sum() == Approx(1.0));

        REQUIRE(dN.col(0).sum() == Approx(0.0).margin(ZERO_MARGIN));
        REQUIRE(dN.col(1).sum() == Approx(0.0).margin(ZERO_MARGIN));
        REQUIRE(dN.col(2).sum() == Approx(0.0).margin(ZERO_MARGIN));
    }
}

TEST_CASE("tetrahedron quadrature")
{
    using namespace neon::quadrature::tetrahedron;

    SECTION("jinyun")
    {
        for (auto&& degree : {1, 2, 3, 4, 5, 6})
        {
            jinyun t{degree};
            REQUIRE(t.degree() >= degree);
            REQUIRE(std::accumulate(begin(t.weights()), end(t.weights()), 0.0) == Approx(1.0 / 6.0));
        }
        REQUIRE(jinyun{1}.points() == 1);
        REQUIRE(jinyun{2}.points() == 4);
        REQUIRE(jinyun{3}.points() == 5);
        REQUIRE(jinyun{4}.points() == 16);
        REQUIRE(jinyun{5}.points() == 17);
        REQUIRE(jinyun{6}.points() == 29);
    }
    SECTION("witherden_vincent")
    {
        for (auto&& degree : {1, 2, 3, 4, 5})
        {
            witherden_vincent t{degree};
            REQUIRE(t.degree() >= degree);
            REQUIRE(std::accumulate(begin(t.weights()), end(t.weights()), 0.0) == Approx(1.0 / 6.0));
        }
        REQUIRE(witherden_vincent{1}.points() == 1);
        REQUIRE(witherden_vincent{2}.points() == 4);
        REQUIRE(witherden_vincent{3}.points() == 8);
        REQUIRE(witherden_vincent{4}.points() == 14);
        REQUIRE(witherden_vincent{5}.points() == 14);
        REQUIRE(witherden_vincent{6}.points() == 24);
    }
}
TEST_CASE("tetrahedron shape functions")
{
    using namespace neon::quadrature::tetrahedron;

    REQUIRE(tetrahedron4{}.number_of_nodes() == 4);
    REQUIRE(tetrahedron10{}.number_of_nodes() == 10);

    SECTION("tetrahedron jinyun quadrature", "[quadrature.tetrahedron.jinyun]")
    {
        for (auto&& degree : {1, 2, 3, 4, 5})
        {
            check_shape_functions(tetrahedron4{}, jinyun{degree});
            check_shape_functions(tetrahedron10{}, jinyun{degree});
        }
    }
    SECTION("tetrahedron witherden vincent", "[quadrature.tetrahedron.witherden_vincent]")
    {
        for (auto&& degree : {1, 2, 3, 4, 5})
        {
            check_shape_functions(tetrahedron4{}, witherden_vincent{degree});
            check_shape_functions(tetrahedron10{}, witherden_vincent{degree});
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
