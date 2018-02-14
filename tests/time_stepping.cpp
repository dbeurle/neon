
#define CATCH_CONFIG_MAIN

#include "catch.hpp"

#include "solver/adaptive_time_step.hpp"
#include "solver/time_step_control.hpp"

#include "io/json.hpp"

using namespace neon;

TEST_CASE("adaptive time control")
{
    json time_data = {{"Period", 1.0},
                      {"Increments", {{"Initial", 1.0}, {"Minimum", 0.1}, {"Maximum", 1.0}}}};

    adaptive_time_step load(time_data, {0.0, 1.0});

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
TEST_CASE("Simple time control")
{
    SECTION("exceptional behaviour")
    {
        REQUIRE_THROWS_AS(time_step_control({{"Ed", 1.0}, {"Start", 0.0}, {"StepSize", 0.1}}),
                          std::domain_error);
        REQUIRE_THROWS_AS(time_step_control({{"End", 1.0}, {"Stat", 0.0}, {"StepSize", 0.1}}),
                          std::domain_error);
        REQUIRE_THROWS_AS(time_step_control({{"End", 1.0}, {"Start", 0.0}, {"Stepize", 0.1}}),
                          std::domain_error);
    }
    SECTION("increments")
    {
        time_step_control times({{"End", 1.0}, {"Start", 0.0}, {"StepSize", 0.1}});

        REQUIRE(times.current_time_step_size() == Approx(0.1));
        REQUIRE(times.number_of_time_steps() == 10);
        REQUIRE(times.is_finished() == false);

        for (auto i = 0; i < 10; i++)
        {
            times.increment();
        }
        REQUIRE(times.is_finished());
    }
}
