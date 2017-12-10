
#define CATCH_CONFIG_MAIN

#include "catch.hpp"

#include "solver/AdaptiveLoadStep.hpp"

#include <string>

#include <json/json.h>

Json::CharReaderBuilder reader;
JSONCPP_STRING input_errors;

using namespace neon;

std::string json_input_file()
{
    return "{\"Name\": \"steel\", \"ElasticModulus\": 1.0, \"PoissonsRatio\": 2.0}";
}

std::string increment_json()
{
    return "{\"Period\" : 1.0, \"Increments\": { "
           "\"Initial\" : 1.0, \"Minimum\" : 0.1, \"Maximum\" : 1.0 }}";
}

TEST_CASE("AdaptiveIncrement")
{
    Json::Value time_data;

    std::istringstream time_data_stream(increment_json());

    REQUIRE(Json::parseFromStream(reader, time_data_stream, &time_data, &input_errors));

    AdaptiveLoadStep load(time_data, {0.0, 1.0});

    SECTION("Basic operation")
    {
        REQUIRE(!load.is_fully_applied());
        REQUIRE(load.step_time() == Approx(1.0));
    }
    SECTION("Well behaved nonlinear iteration")
    {
        load.update_convergence_state(true);
        REQUIRE(load.is_fully_applied());
    }
    SECTION("Badly behaved nonlinear iteration")
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
