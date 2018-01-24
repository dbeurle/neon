
#define CATCH_CONFIG_MAIN

#include "catch.hpp"

#include "solver/AdaptiveLoadStep.hpp"

#include "io/json.hpp"

using namespace neon;

TEST_CASE("adaptive increment")
{
    json time_data = {{"Period", 1.0},
                      {"Increments", {{"Initial", 1.0}, {"Minimum", 0.1}, {"Maximum", 1.0}}}};

    AdaptiveLoadStep load(time_data, {0.0, 1.0});

    SECTION("basic operation")
    {
        REQUIRE(!load.is_fully_applied());
        REQUIRE(load.step_time() == Approx(1.0));
    }
    SECTION("well behaved nonlinear iteration")
    {
        load.update_convergence_state(true);
        REQUIRE(load.is_fully_applied());
    }
    SECTION("naughty nonlinear iteration")
    {
        auto last_good_load_factor = load.step_time();
        for (int i = 0; i < 9; i++)
        {
            load.update_convergence_state(false);
            REQUIRE(load.step_time() < last_good_load_factor);
        }
        REQUIRE_THROWS_AS(load.update_convergence_state(false), std::runtime_error);
    }
}
