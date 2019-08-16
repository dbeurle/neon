
#include <catch2/catch.hpp>

#include "quadrature/hexahedron/gauss_legendre.hpp"
#include "interpolations/hexahedron.hpp"

#include "mesh/element_topology.hpp"
#include "interpolations/interpolation_factory.hpp"

#include <numeric>

constexpr auto ZERO_MARGIN = 1.0e-5;

using namespace neon;

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

TEST_CASE("Hexahedron quadrature scheme test", "[quadrature.hexahedron]")
{
    using namespace neon::quadrature::hexahedron;

    SECTION("Hexahedron", "[quadrature.hexahedron.gauss_legendre]")
    {
        for (auto&& degree : {1, 2, 3, 4, 5})
        {
            gauss_legendre h(degree);
            REQUIRE(h.degree() >= degree);
            REQUIRE(std::accumulate(begin(h.weights()), end(h.weights()), 0.0) == Approx(8.0));
        }
        REQUIRE(gauss_legendre{1}.points() == 1);
        REQUIRE(gauss_legendre{2}.points() == 8);
        REQUIRE(gauss_legendre{3}.points() == 8);
        REQUIRE(gauss_legendre{4}.points() == 27);
        REQUIRE(gauss_legendre{5}.points() == 27);
    }
    SECTION("hexahedron shape functions", "[quadrature.hexahedron.gauss_legendre]")
    {
        REQUIRE(hexahedron8{}.number_of_nodes() == 8);
        REQUIRE(hexahedron20{}.number_of_nodes() == 20);
        REQUIRE(hexahedron27{}.number_of_nodes() == 27);

        for (auto&& degree : {1, 2, 3, 4, 5})
        {
            check_shape_functions(hexahedron8{}, gauss_legendre{degree});
            check_shape_functions(hexahedron20{}, gauss_legendre{degree});
            check_shape_functions(hexahedron27{}, gauss_legendre{degree});
        }
    }
}
// TEST_CASE("Hexahedron volume evaluation")
// {
//     SECTION("hexahedron8 volume evaluation")
//     {
//         hexahedron8 hex8();
//
//         matrix x(3, 8);
//         x << 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, //
//             0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0,  //
//             0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0;
//
//         REQUIRE(hex8.compute_measure(x) == Approx(1.0));
//     }
//     SECTION("hexahedron20 volume evaluation")
//     {
//         hexahedron20 hex20;
//
//         // xyz coordinates of the unit cube
//         matrix x(3, 20);
//         x << 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.5, 1.0, 0.5, 0.0, 0.5, 1.0, 0.5, 0.0, 0.0,
//             1.0, 1.0, 0.0, //
//             0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, /**/ 0.0, 0.5, 1.0, 0.5, 0.0, 0.5, 1.0, 0.5,
//             0.0, 0.0, 1.0, 1.0, //
//             0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 0.5,
//             0.5, 0.5, 0.5;
//
//         REQUIRE(hex20.compute_measure(x) == Approx(1.0));
//     }
//     SECTION("hexahedron27 volume evaluation")
//     {
//         SECTION("Six point rule")
//         {
//             hexahedron27 hex27;
//
//             // xyz coordinates of the unit cube
//             matrix x(3, 27);
//             x << 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, /**/ 0.5, 1.0, 0.5, 0.0, 0.5, 1.0, 0.5,
//                 0.0, 0.0, 1.0, 1.0, 0.0, /**/ 0.5, 0.5, 0.5, 0.5, 0.0, 1.0, 0.5, //
//                 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, /**/ 0.0, 0.5, 1.0, 0.5, 0.0, 0.5, 1.0, 0.5,
//                 0.0, 0.0, 1.0, 1.0, /**/ 0.5, 0.5, 0.0, 1.0, 0.5, 0.5, 0.5, //
//                 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, /**/ 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0,
//                 0.5, 0.5, 0.5, 0.5, /**/ 0.0, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5;
//
//             REQUIRE(hex27.compute_measure(x) == Approx(1.0));
//
//             x *= 2.0;
//
//             REQUIRE(hex27.compute_measure(x) == Approx(8.0));
//         }
//         SECTION("Eight point rule")
//         {
//             hexahedron27 hex27;
//
//             // xyz coordinates of the unit cube
//             matrix x(3, 27);
//             x << 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, /**/ 0.5, 1.0, 0.5, 0.0, 0.5, 1.0, 0.5,
//                 0.0, 0.0, 1.0, 1.0, 0.0, /**/ 0.5, 0.5, 0.5, 0.5, 0.0, 1.0, 0.5, //
//                 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, /**/ 0.0, 0.5, 1.0, 0.5, 0.0, 0.5, 1.0, 0.5,
//                 0.0, 0.0, 1.0, 1.0, /**/ 0.5, 0.5, 0.0, 1.0, 0.5, 0.5, 0.5, //
//                 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, /**/ 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0,
//                 0.5, 0.5, 0.5, 0.5, /**/ 0.0, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5;
//
//             REQUIRE(hex27.compute_measure(x) == Approx(1.0));
//
//             x *= 2.0;
//
//             REQUIRE(hex27.compute_measure(x) == Approx(8.0));
//         }
//     }
// }
