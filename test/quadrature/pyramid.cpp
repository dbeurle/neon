
#include <catch2/catch.hpp>

#include "interpolations/pyramid.hpp"
#include "quadrature/pyramid/bedrosian.hpp"

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

// static matrix3x unit_pyramid5_coordinates()
// {
//     auto constexpr base = 1.0;
//     auto constexpr height = 1.0;
//
//     matrix3x x(3, 5);
//     x << 0.0, base, base, 0.0, height / 2.0, //
//         0.0, 0.0, base, base, height / 2.0,  //
//         0.0, 0.0, 0.0, 0.0, height;
//     return x;
// }

// static matrix3x unit_pyramid13_coordinates()
// {
//     matrix3x x(3, 13);
//     x << 1.0, 0.0, -1.0, -1.0, -1.0, 0.0, 1.0, 1.0, 0.5, -0.5, -0.5, 0.5, 0.0, //
//         1.0, 1.0, 1.0, 0.0, -1.0, -1.0, -1.0, 0.0, 0.5, 0.5, -0.5, -0.5, 0.0,  //
//         0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.5, 0.5, 0.5, 0.5, 1.0;
//
//     REQUIRE(element.compute_measure(x) == Approx(4.0 / 3.0));
// }

TEST_CASE("pyramid")
{
    using namespace neon::quadrature::pyramid;

    SECTION("bedrosian")
    {
        for (auto&& degree : {1, 2, 3, 4, 5})
        {
            bedrosian p{degree};
            REQUIRE(p.degree() >= degree);
            REQUIRE(std::accumulate(begin(p.weights()), end(p.weights()), 0.0) == Approx(4.0 / 3.0));
        }
        REQUIRE(bedrosian{1}.points() == 1);
        REQUIRE(bedrosian{2}.points() == 8);
        REQUIRE(bedrosian{5}.points() == 27);
    }
    SECTION("prism shape functions")
    {
        pyramid5 element;
        bedrosian scheme(1);

        REQUIRE(pyramid5{}.number_of_nodes() == 5);
        REQUIRE(pyramid13{}.number_of_nodes() == 13);

        for (auto&& degree : {1, 2, 3, 4, 5})
        {
            check_shape_functions(pyramid5{}, bedrosian{degree});
            check_shape_functions(pyramid13{}, bedrosian{degree});
        }
    }
    // SECTION("Five node volume - one point")
    // {
    //     pyramid5 element(pyramid_quadrature::point::one);
    //
    //     auto const base = 1.0;
    //     auto const height = 1.0;
    //
    //
    //     REQUIRE(element.compute_measure(x) == Approx(1.0 / 3.0));
    // }
    // SECTION("Five node volume - eight point")
    // {
    //     pyramid5 element(pyramid_quadrature::point::eight);
    //
    //     auto constexpr base = 1.0;
    //     auto constexpr height = 1.0;
    //
    //     matrix3x x(3, 5);
    //     x << 0.0, base, base, 0.0, height / 2.0, //
    //         0.0, 0.0, base, base, height / 2.0,  //
    //         0.0, 0.0, 0.0, 0.0, height;
    //
    //     REQUIRE(element.compute_measure(x) == Approx(1.0 / 3.0));
    // }
    // SECTION("Thirteen node volume - eight point")
    // {
    //     pyramid13 element(pyramid_quadrature::point::eight);
    //
    //     matrix3x x(3, 13);
    //     x << 1.0, 0.0, -1.0, -1.0, -1.0, 0.0, 1.0, 1.0, 0.5, -0.5, -0.5, 0.5, 0.0, //
    //         1.0, 1.0, 1.0, 0.0, -1.0, -1.0, -1.0, 0.0, 0.5, 0.5, -0.5, -0.5, 0.0,  //
    //         0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.5, 0.5, 0.5, 0.5, 1.0;
    //
    //     REQUIRE(element.compute_measure(x) == Approx(4.0 / 3.0));
    // }
}
