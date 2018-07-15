
#include <catch.hpp>

#include "geometry/profile_factory.hpp"
#include "io/json.hpp"

#include <stdexcept>

using namespace neon;

TEST_CASE("Geometry factory")
{
    SECTION("Name validation")
    {
        REQUIRE_THROWS_AS(geometry::make_profile(json::parse("{\"proile\": \"rectangle\"}")),
                          std::domain_error);
        REQUIRE_THROWS_AS(geometry::make_profile(json::parse("{\"profile\": \"rect\"}")),
                          std::domain_error);
    }
}
TEST_CASE("Dimensional values")
{
    SECTION("Rectangle")
    {
        REQUIRE_THROWS_AS(geometry::make_profile(
                              json{{"profile", "rectangle"}, {"wdth", 2.0}, {"height", 0.3}}),
                          std::domain_error);

        REQUIRE_THROWS_AS(geometry::make_profile(
                              json{{"profile", "rectangle"}, {"width", 2.0}, {"hight", 0.3}}),
                          std::domain_error);

        REQUIRE_THROWS_AS(geometry::make_profile(
                              json{{"profile", "rectangle"}, {"width", -2.0}, {"height", 0.3}}),
                          std::domain_error);

        REQUIRE_THROWS_AS(geometry::make_profile(
                              json{{"profile", "rectangle"}, {"width", 2.0}, {"height", -0.3}}),
                          std::domain_error);
    }
    SECTION("Circle")
    {
        REQUIRE_THROWS_AS(geometry::make_profile(json{{"profile", "circle"}, {"dimeter", 2.0}}),
                          std::domain_error);

        REQUIRE_THROWS_AS(geometry::make_profile(json{{"profile", "circle"}, {"diameter", -2.0}}),
                          std::domain_error);
    }
}

TEST_CASE("Calculated values")
{
    SECTION("Square")
    {
        auto profile = geometry::make_profile(
            json{{"profile", "rectangle"}, {"width", 1.0}, {"height", 1.0}});

        REQUIRE(profile->area() == Approx(1.0));
        REQUIRE(profile->second_moment_area().first == Approx(1.0 / 12.0));
        REQUIRE(profile->second_moment_area().second == Approx(1.0 / 12.0));
        REQUIRE(profile->shear_area().first == Approx(1.0));
        REQUIRE(profile->shear_area().second == Approx(1.0));
    }
    SECTION("Horizontal rectangle")
    {
        auto profile = geometry::make_profile(
            json{{"profile", "rectangle"}, {"width", 2.0}, {"height", 1.0}});

        REQUIRE(profile->area() == Approx(2.0));
        REQUIRE(profile->second_moment_area().first == Approx(1.0 / 6.0));
        REQUIRE(profile->second_moment_area().second == Approx(2.0 / 3.0));
        REQUIRE(profile->shear_area().first == Approx(2.0));
        REQUIRE(profile->shear_area().second == Approx(2.0));
    }
    SECTION("Vertical rectangle")
    {
        auto profile = geometry::make_profile(
            json{{"profile", "rectangle"}, {"width", 1.0}, {"height", 2.0}});

        REQUIRE(profile->area() == Approx(2.0));
        REQUIRE(profile->second_moment_area().first == Approx(2.0 / 3.0));
        REQUIRE(profile->second_moment_area().second == Approx(1.0 / 6.0));
        REQUIRE(profile->shear_area().first == Approx(2.0));
        REQUIRE(profile->shear_area().second == Approx(2.0));
    }
    SECTION("Small circle")
    {
        auto profile = geometry::make_profile(json{{"profile", "circle"}, {"diameter", 1.0}});

        REQUIRE(profile->area() == Approx(0.785398));
        REQUIRE(profile->second_moment_area().first == Approx(0.049087));
        REQUIRE(profile->second_moment_area().second == Approx(0.049087));
        REQUIRE(profile->shear_area().first == Approx(0.785398));
        REQUIRE(profile->shear_area().second == Approx(0.785398));
    }
    SECTION("Large circle")
    {
        auto profile = geometry::make_profile(json{{"profile", "circle"}, {"diameter", 2.0}});

        REQUIRE(profile->area() == Approx(3.141593));
        REQUIRE(profile->second_moment_area().first == Approx(0.785398));
        REQUIRE(profile->second_moment_area().second == Approx(0.785398));
        REQUIRE(profile->shear_area().first == Approx(3.141593));
        REQUIRE(profile->shear_area().second == Approx(3.141593));
    }
}
