#define CATCH_CONFIG_MAIN // This tells Catch to provide a main() - only do this in one
                          // cpp file
#include "catch.hpp"

#include "solver/AdaptiveLoadStep.hpp"

#include <string>

#include <json/json.h>

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
    Json::Reader time_file;

    REQUIRE(time_file.parse(increment_json().c_str(), time_data));

    AdaptiveLoadStep load(time_data);

    SECTION("Basic operation")
    {
        REQUIRE(!load.is_fully_applied());
        REQUIRE(load.factor() == Approx(1.0));
    }
    SECTION("Well behaved nonlinear iteration")
    {
        load.update_convergence_state(true);
        REQUIRE(load.is_fully_applied());
    }
    SECTION("Badly behaved nonlinear iteration")
    {
        auto last_good_load_factor = load.factor();
        for (int i = 0; i < 9; i++)
        {
            load.update_convergence_state(false);
            REQUIRE(load.factor() < last_good_load_factor);
        }
        REQUIRE_THROWS_AS(load.update_convergence_state(false), std::runtime_error);
    }
}
