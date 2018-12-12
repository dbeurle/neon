
#include <catch2/catch.hpp>

#include "quadrature/sphere/unit_sphere_quadrature.hpp"

#include <range/v3/numeric/accumulate.hpp>
#include <numeric>

using namespace neon;

constexpr auto ZERO_MARGIN = 1.0e-5;

TEST_CASE("Unit sphere quadrature scheme test", "[unit_sphere_quadrature]")
{
    SECTION("BO21 rule test")
    {
        unit_sphere_quadrature unit_sphere(unit_sphere_quadrature::point::BO21);

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
        unit_sphere_quadrature unit_sphere(unit_sphere_quadrature::point::BO33);

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
        unit_sphere_quadrature unit_sphere(unit_sphere_quadrature::point::FM900);

        REQUIRE(unit_sphere.points() == 900);
        REQUIRE(std::accumulate(begin(unit_sphere.weights()), end(unit_sphere.weights()), 0.0)
                == Approx(1.0));

        for (auto const& coordinate : unit_sphere.coordinates())
        {
            auto const& [l, x, y, z] = coordinate;
            vector3 norm_check(x, y, z);
            REQUIRE(norm_check.norm() == Approx(1.0));
        }
    }
}
