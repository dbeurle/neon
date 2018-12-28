
#include <catch2/catch.hpp>

#include "geometry/profile_factory.hpp"
#include "io/json.hpp"

#include <stdexcept>

using namespace neon;

TEST_CASE("Geometry factory")
{
    SECTION("Name validation")
    {
        REQUIRE_THROWS_AS(geometry::make_profile(json::parse("{\"tpe\": \"rectangle\"}")),
                          std::domain_error);
        REQUIRE_THROWS_AS(geometry::make_profile(json::parse("{\"type\": \"rect\"}")),
                          std::domain_error);
    }
}
TEST_CASE("Dimensional values")
{
    SECTION("Rectangle")
    {
        // Spelling
        REQUIRE_THROWS_AS(geometry::make_profile(
                              json{{"type", "rectangle"}, {"wdth", 2.0}, {"height", 0.3}}),
                          std::domain_error);

        REQUIRE_THROWS_AS(geometry::make_profile(
                              json{{"type", "rectangle"}, {"width", 2.0}, {"hight", 0.3}}),
                          std::domain_error);

        // Negative value
        REQUIRE_THROWS_AS(geometry::make_profile(
                              json{{"type", "rectangle"}, {"width", -2.0}, {"height", 0.3}}),
                          std::domain_error);

        REQUIRE_THROWS_AS(geometry::make_profile(
                              json{{"type", "rectangle"}, {"width", 2.0}, {"height", -0.3}}),
                          std::domain_error);
    }
    SECTION("Circle")
    {
        REQUIRE_THROWS_AS(geometry::make_profile(json{{"type", "circle"}, {"dimeter", 2.0}}),
                          std::domain_error);

        REQUIRE_THROWS_AS(geometry::make_profile(json{{"type", "circle"}, {"diameter", -2.0}}),
                          std::domain_error);
    }
    SECTION("Hollow circle")
    {
        // Spelling
        REQUIRE_THROWS_AS(geometry::make_profile(json{{"type", "hollow_circle"},
                                                      {"iner_diameter", 1.0},
                                                      {"outer_diameter", 2.0}}),
                          std::domain_error);
        REQUIRE_THROWS_AS(geometry::make_profile(json{{"type", "hollow_circle"},
                                                      {"inner_diameter", 1.0},
                                                      {"outer_dimeter", 2.0}}),
                          std::domain_error);
        // Negative Value
        REQUIRE_THROWS_AS(geometry::make_profile(json{{"type", "hollow_circle"},
                                                      {"inner_diameter", -1.0},
                                                      {"outer_diameter", 2.0}}),
                          std::domain_error);
        REQUIRE_THROWS_AS(geometry::make_profile(json{{"type", "hollow_circle"},
                                                      {"inner_diameter", 1.0},
                                                      {"outer_diameter", -2.0}}),
                          std::domain_error);
        // Invalid inner_diameter and outer_diameter inputs
        REQUIRE_THROWS_AS(geometry::make_profile(json{{"type", "hollow_circle"},
                                                      {"inner_diameter", 2.0},
                                                      {"outer_diameter", 2.0}}),
                          std::domain_error);
        REQUIRE_THROWS_AS(geometry::make_profile(json{{"type", "hollow_circle"},
                                                      {"inner_diameter", 2.0},
                                                      {"outer_diameter", 1.0}}),
                          std::domain_error);
    }
    SECTION("Rectangular angle")
    {
        // Spelling
        REQUIRE_THROWS_AS(geometry::make_profile(json{{"type", "rectangular_angle"},
                                                      {"with", 1.0},
                                                      {"height", 2.0},
                                                      {"thickness", 0.5}}),
                          std::domain_error);
        REQUIRE_THROWS_AS(geometry::make_profile(json{{"type", "rectangular_angle"},
                                                      {"width", 1.0},
                                                      {"heigt", 2.0},
                                                      {"thickness", 0.5}}),
                          std::domain_error);
        REQUIRE_THROWS_AS(geometry::make_profile(json{{"type", "rectangular_angle"},
                                                      {"width", 1.0},
                                                      {"height", 2.0},
                                                      {"thickess", 0.5}}),
                          std::domain_error);
        // Negative value
        REQUIRE_THROWS_AS(geometry::make_profile(json{{"type", "rectangular_angle"},
                                                      {"width", -1.0},
                                                      {"height", 2.0},
                                                      {"thickess", 0.5}}),
                          std::domain_error);
        REQUIRE_THROWS_AS(geometry::make_profile(json{{"type", "rectangular_angle"},
                                                      {"width", 1.0},
                                                      {"height", -2.0},
                                                      {"thickess", 0.5}}),
                          std::domain_error);
        REQUIRE_THROWS_AS(geometry::make_profile(json{{"type", "rectangular_angle"},
                                                      {"width", 1.0},
                                                      {"height", 2.0},
                                                      {"thickess", -0.5}}),
                          std::domain_error);
    }
}

TEST_CASE("Calculated values")
{
    SECTION("Square")
    {
        auto profile = geometry::make_profile(
            json{{"type", "rectangle"}, {"width", 1.0}, {"height", 1.0}});

        REQUIRE(profile->area() == Approx(1.0));
        REQUIRE(profile->second_moment_area().first == Approx(1.0 / 12.0));
        REQUIRE(profile->second_moment_area().second == Approx(1.0 / 12.0));
        REQUIRE(profile->shear_area().first == Approx(1.0));
        REQUIRE(profile->shear_area().second == Approx(1.0));
    }
    SECTION("Horizontal rectangle")
    {
        auto profile = geometry::make_profile(
            json{{"type", "rectangle"}, {"width", 2.0}, {"height", 1.0}});

        REQUIRE(profile->area() == Approx(2.0));
        REQUIRE(profile->second_moment_area().first == Approx(1.0 / 6.0));
        REQUIRE(profile->second_moment_area().second == Approx(2.0 / 3.0));
        REQUIRE(profile->shear_area().first == Approx(2.0));
        REQUIRE(profile->shear_area().second == Approx(2.0));
    }
    SECTION("Vertical rectangle")
    {
        auto profile = geometry::make_profile(
            json{{"type", "rectangle"}, {"width", 1.0}, {"height", 2.0}});

        REQUIRE(profile->area() == Approx(2.0));
        REQUIRE(profile->shear_area().first == Approx(2.0));
        REQUIRE(profile->shear_area().second == Approx(2.0));
        REQUIRE(profile->second_moment_area().first == Approx(2.0 / 3.0));
        REQUIRE(profile->second_moment_area().second == Approx(1.0 / 6.0));
    }
    SECTION("Small circle")
    {
        auto profile = geometry::make_profile(json{{"type", "circle"}, {"diameter", 1.0}});

        REQUIRE(profile->area() == Approx(0.785398));
        REQUIRE(profile->shear_area().first == Approx(0.785398));
        REQUIRE(profile->shear_area().second == Approx(0.785398));
        REQUIRE(profile->second_moment_area().first == Approx(0.049087));
        REQUIRE(profile->second_moment_area().second == Approx(0.049087));
    }
    SECTION("Large circle")
    {
        auto profile = geometry::make_profile(json{{"type", "circle"}, {"diameter", 2.0}});

        REQUIRE(profile->area() == Approx(3.141593));
        REQUIRE(profile->shear_area().first == Approx(3.141593));
        REQUIRE(profile->shear_area().second == Approx(3.141593));
        REQUIRE(profile->second_moment_area().first == Approx(0.785398));
        REQUIRE(profile->second_moment_area().second == Approx(0.785398));
    }
    SECTION("Hollow circle")
    {
        auto profile = geometry::make_profile(
            json{{"type", "hollow_circle"}, {"inner_diameter", 1.0}, {"outer_diameter", 2.0}});

        REQUIRE(profile->area() == Approx(2.35619));
        REQUIRE(profile->shear_area().first == Approx(2.35619));
        REQUIRE(profile->shear_area().second == Approx(2.35619));
        REQUIRE(profile->second_moment_area().first == Approx(0.736311));
        REQUIRE(profile->second_moment_area().second == Approx(0.736311));
    }
    SECTION("Rectangular angle")
    {
        auto profile = geometry::make_profile(
            json{{"type", "rectangular_angle"}, {"width", 1.0}, {"height", 2.0}, {"thickness", 0.5}});

        REQUIRE(profile->area() == Approx(1.25));
        REQUIRE(profile->second_moment_area().first == Approx(0.451042));
        REQUIRE(profile->second_moment_area().second == Approx(0.076042));
    }
    {
        auto profile = geometry::make_profile(
            json{{"type", "rectangular_angle"}, {"width", 2.0}, {"height", 1.0}, {"thickness", 0.5}});

        REQUIRE(profile->area() == Approx(1.25));
        REQUIRE(profile->second_moment_area().first == Approx(0.076042));
        REQUIRE(profile->second_moment_area().second == Approx(0.451042));
    }
}
