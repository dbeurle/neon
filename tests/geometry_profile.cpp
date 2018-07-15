
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
