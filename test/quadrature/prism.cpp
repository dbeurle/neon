
#include <catch2/catch.hpp>

#include "quadrature/prism/felippa.hpp"
#include "quadrature/prism/kubatko.hpp"

#include "interpolations/prism.hpp"

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

static matrix3x six_node_coordinates()
{
    matrix3x x(3, 6);
    x << 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, //
        0.0, 1.0, 0.0, 0.0, 1.0, 0.0,  //
        -1.0, -1.0, -1.0, 1.0, 1.0, 1.0;
    return x;
}

static matrix3x fifteen_node_coordinates()
{
    matrix3x x(3, 15);
    x << 1.0, 0.0, 0.0, 0.5, 0.0, 0.5, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.5, 0.0, 0.5, //
        0.0, 1.0, 0.0, 0.5, 0.5, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.5, 0.5, 0.0,  //
        -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0;
    return x;
}

TEST_CASE("Prism quadrature scheme test")
{
    using namespace neon::quadrature::prism;

    SECTION("Felippa quadrature")
    {
        for (auto&& degree : {1, 2, 3, 4, 5, 6})
        {
            felippa p{degree};
            REQUIRE(p.degree() >= degree);
            REQUIRE(std::accumulate(begin(p.weights()), end(p.weights()), 0.0) == Approx(1.0));
        }
        REQUIRE(felippa{1}.points() == 1);
        REQUIRE(felippa{2}.points() == 6);
        REQUIRE(felippa{3}.points() == 6);
        REQUIRE(felippa{4}.points() == 18);
        REQUIRE(felippa{5}.points() == 21);
        REQUIRE(felippa{6}.points() == 48);
    }
    SECTION("Kubatko quadrature")
    {
        for (auto&& degree : {1, 2, 3, 4, 5})
        {
            kubatko p{degree};
            REQUIRE(p.degree() >= degree);
            REQUIRE(std::accumulate(begin(p.weights()), end(p.weights()), 0.0) == Approx(1.0));
        }
        REQUIRE(kubatko{1}.points() == 1);
        REQUIRE(kubatko{2}.points() == 4);
        REQUIRE(kubatko{3}.points() == 6);
        REQUIRE(kubatko{4}.points() == 11);
        REQUIRE(kubatko{5}.points() == 15);
    }
}
TEST_CASE("Element checks")
{
    using namespace neon::quadrature::prism;

    REQUIRE(prism6{}.number_of_nodes() == 6);
    REQUIRE(prism15{}.number_of_nodes() == 15);

    SECTION("prism shape functions")
    {
        for (auto&& degree : {1, 2, 3, 4, 5})
        {
            check_shape_functions(prism6{}, kubatko{degree});
            check_shape_functions(prism6{}, felippa{degree});
            check_shape_functions(prism15{}, kubatko{degree});
            check_shape_functions(prism15{}, felippa{degree});
        }
    }
    // SECTION("Six node - six point area")
    // {
    //     prism6 pri6(prism_quadrature::point::six);
    //
    //     REQUIRE(pri6.compute_measure(six_node_coordinates()) == Approx(1.0));
    // }
    // SECTION("Fifteen node - six point evaluation")
    // {
    //     prism15 element(prism_quadrature::point::six);
    //
    //     REQUIRE(element.number_of_nodes() == 15);
    //     REQUIRE(element.quadrature().points() == 6);
    //
    //     element.quadrature().for_each([&](auto const& femval, auto) {
    //         auto const [N, dN] = femval;
    //
    //         REQUIRE(N.sum() == Approx(1.0));
    //
    //         REQUIRE(dN.col(0).sum() == Approx(0.0).margin(ZERO_MARGIN));
    //         REQUIRE(dN.col(1).sum() == Approx(0.0).margin(ZERO_MARGIN));
    //         REQUIRE(dN.col(2).sum() == Approx(0.0).margin(ZERO_MARGIN));
    //     });
    // }
    // SECTION("Fifteen node - nine point area")
    // {
    //     prism15 pri15(prism_quadrature::point::nine);
    //
    //     REQUIRE(pri15.compute_measure(fifteen_node_coordinates()) == Approx(1.0));
    // }
}
