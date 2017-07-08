#define CATCH_CONFIG_MAIN // This tells Catch to provide a main() - only do this in one cpp file
#include "catch.hpp"

#include "solver/AdaptiveIncrement.hpp"

#include <json/json.h>

using namespace neon;

std::string json_input_file()
{
    return "{\"Name\": \"steel\", \"ElasticModulus\": 1.0, \"PoissonsRatio\": 2.0}";
}

TEST_CASE("AdaptiveIncrement")
{
    AdaptiveIncrement load(0.1, 1.0, 100);

    REQUIRE(!load.is_fully_applied());
    REQUIRE(load.factor() == Approx(0.1));

    double load_factor = 0.1;

    // Simulate a well behaved nonlinear-iteration
    SECTION("Well behaved nonlinear iteration")
    {
        for (int i = 0; i < 3; i++)
        {
            load.update_convergence_state(true);

            REQUIRE(!load.is_fully_applied());
            REQUIRE(load.factor() > load_factor);
            REQUIRE(load.factor() <= Approx(1.0));

            load_factor = load.factor();
        }
        load.update_convergence_state(true);

        REQUIRE(load.is_fully_applied());
    }
    SECTION("Normally behaved nonlinear iteration")
    {
        auto last_good_load_factor = load.factor();
        load.update_convergence_state(true);

        for (int i = 0; i < 4; i++)
        {
            load.update_convergence_state(false);
            REQUIRE(load.factor() > last_good_load_factor);
        }
    }
    SECTION("Badly behaved nonlinear iteration")
    {
        auto last_good_load_factor = load.factor();
        load.update_convergence_state(true);

        for (int i = 0; i < 4; i++)
        {
            load.update_convergence_state(false);
            REQUIRE(load.factor() > last_good_load_factor);
        }
        REQUIRE_THROWS_AS(load.update_convergence_state(false), std::runtime_error);
    }
}
