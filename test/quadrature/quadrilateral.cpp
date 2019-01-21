
#include <catch2/catch.hpp>

#include "quadrature/quadrilateral/gauss_legendre.hpp"
#include "interpolations/quadrilateral.hpp"

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
    }
}

TEST_CASE("Quadrilateral quadrature scheme test")
{
    using namespace neon::quadrature::quadrilateral;

    SECTION("quadrilateral gauss legendre", "[quadrature.quadrilateral.gauss_legendre]")
    {
        for (auto&& degree : {1, 2, 3, 4, 5})
        {
            gauss_legendre h(degree);
            REQUIRE(h.degree() >= degree);
            REQUIRE(std::accumulate(begin(h.weights()), end(h.weights()), 0.0) == Approx(4.0));
        }
        REQUIRE(gauss_legendre{1}.points() == 1);
        REQUIRE(gauss_legendre{2}.points() == 4);
        REQUIRE(gauss_legendre{3}.points() == 4);
        REQUIRE(gauss_legendre{4}.points() == 9);
        REQUIRE(gauss_legendre{5}.points() == 9);
    }
    SECTION("quadrilateral shape functions", "[quadrature.quadrilateral.gauss_legendre]")
    {
        REQUIRE(quadrilateral4{}.number_of_nodes() == 4);
        REQUIRE(quadrilateral8{}.number_of_nodes() == 8);
        REQUIRE(quadrilateral9{}.number_of_nodes() == 9);

        for (auto&& degree : {1, 2, 3, 4, 5})
        {
            check_shape_functions(quadrilateral4{}, gauss_legendre{degree});
            check_shape_functions(quadrilateral8{}, gauss_legendre{degree});
            check_shape_functions(quadrilateral9{}, gauss_legendre{degree});
        }
    }
    // SECTION("quadrilateral4 surface area - one point")
    // {
    //     quadrilateral4 quad4;
    //
    //     matrix x(3, 4);
    //     x << 0.0, 1.0, 1.0, 0.0, //
    //         0.0, 0.0, 1.0, 1.0,  //
    //         0.0, 0.0, 0.0, 0.0;
    //
    //     REQUIRE(quad4.compute_measure(x) == Approx(1.0));
    // }
    // SECTION("quadrilateral4 surface area - four point")
    // {
    //     quadrilateral4 quad4(quadril);
    //
    //     matrix x(3, 4);
    //     x << 0.0, 1.0, 1.0, 0.0, //
    //         0.0, 0.0, 1.0, 1.0,  //
    //         0.0, 0.0, 0.0, 0.0;
    //
    //     REQUIRE(quad4.compute_measure(x) == Approx(1.0));
    // }
    // SECTION("quadrilateral8 surface area - four point")
    // {
    //     quadrilateral8 quad8(quadril);
    //
    //     matrix x(3, 8);
    //     x << 0.0, 1.0, 1.0, 0.0, 0.5, 1.0, 0.5, 0.0, //
    //         0.0, 0.0, 1.0, 1.0, 0.0, 0.5, 1.0, 0.5,  //
    //         0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0;
    //
    //     REQUIRE(quad8.compute_measure(x) == Approx(1.0));
    // }
    // SECTION("quadrilateral8 interpolation function - nine point")
    // {
    //     quadrilateral8 quad8(quadril);
    //
    //     matrix x(3, 8);
    //     x << 0.0, 1.0, 1.0, 0.0, 0.5, 1.0, 0.5, 0.0, //
    //         0.0, 0.0, 1.0, 1.0, 0.0, 0.5, 1.0, 0.5,  //
    //         0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0;
    //
    //     REQUIRE(quad8.compute_measure(x) == Approx(1.0));
    // }
    // SECTION("quadrilateral9 surface area - four point")
    // {
    //     quadrilateral9 quad9(quadril);
    //
    //     matrix x(3, 9);
    //     x << 0.0, 1.0, 1.0, 0.0, 0.5, 1.0, 0.5, 0.0, 0.5, //
    //         0.0, 0.0, 1.0, 1.0, 0.0, 0.5, 1.0, 0.5, 0.5,  //
    //         0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0;
    //
    //     REQUIRE(quad9.compute_measure(x) == Approx(1.0));
    // }
    // SECTION("quadrilateral9 interpolation function - nine point")
    // {
    //     quadrilateral9 quad9(quadril);
    //
    //     matrix x(3, 9);
    //     x << 0.0, 1.0, 1.0, 0.0, 0.5, 1.0, 0.5, 0.0, 0.5, //
    //         0.0, 0.0, 1.0, 1.0, 0.0, 0.5, 1.0, 0.5, 0.5,  //
    //         0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0;
    //
    //     REQUIRE(quad9.compute_measure(x) == Approx(1.0));
    // }
}
