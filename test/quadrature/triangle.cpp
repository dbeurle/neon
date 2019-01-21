
#include <catch2/catch.hpp>

#include "quadrature/triangle/cowper.hpp"
#include "interpolations/triangle.hpp"

#include "mesh/element_topology.hpp"
#include "interpolations/interpolation_factory.hpp"

#include <numeric>

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
    }
}

using namespace neon;

TEST_CASE("Triangle quadrature")
{
    using namespace neon::quadrature::triangle;

    REQUIRE(triangle3{}.number_of_nodes() == 3);
    REQUIRE(triangle6{}.number_of_nodes() == 6);

    SECTION("triangle")
    {
        for (auto&& degree : {1, 2, 3, 4, 5, 6, 7})
        {
            cowper t(degree);
            REQUIRE(t.degree() >= degree);
            REQUIRE(std::accumulate(begin(t.weights()), end(t.weights()), 0.0) == Approx(1.0 / 2.0));
        }
        REQUIRE(cowper{1}.points() == 3);
        REQUIRE(cowper{2}.points() == 3);
        REQUIRE(cowper{3}.points() == 6);
        REQUIRE(cowper{4}.points() == 6);
        REQUIRE(cowper{5}.points() == 9);
        REQUIRE(cowper{6}.points() == 12);
        REQUIRE(cowper{7}.points() == 13);
        REQUIRE_THROWS_AS(cowper{8}, std::domain_error);
    }
    SECTION("triangle cowper")
    {
        for (auto&& degree : {1, 2, 3, 4, 5, 6, 7})
        {
            check_shape_functions(triangle3{}, cowper{degree});
            check_shape_functions(triangle6{}, cowper{degree});
        }
    }
}
// SECTION("triangle3 surface area")
// {
//     triangle3 tri3(triangle_quadrature::point::one);
//
//     matrix x(3, 3);
//     x << 0.0, 0.0, 0.0, //
//         1.0, 0.0, 0.0,  //
//         0.0, 1.0, 0.0;
//     x.transposeInPlace();
//
//     REQUIRE(tri3.compute_measure(x) == Approx(0.5));
// }
// SECTION("triangle6 surface area")
// {
//     triangle6 patch(triangle_quadrature::point::three);
//
//     matrix x(6, 3);
//     x << 0.0, 0.0, 0.0, //
//         1.0, 0.0, 0.0,  //
//         0.0, 1.0, 0.0,  //
//         0.5, 0.0, 0.0,  //
//         0.5, 0.5, 0.0,  //
//         0.0, 0.5, 0.0;
//     x.transposeInPlace();
//
//     REQUIRE(patch.compute_measure(x) == Approx(0.5));
// }
