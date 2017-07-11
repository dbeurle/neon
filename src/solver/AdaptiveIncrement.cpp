
#include "AdaptiveIncrement.hpp"

#include <exception>
#include <range/v3/algorithm.hpp>
#include <termcolor/termcolor.hpp>

#include <json/value.h>

namespace neon
{
AdaptiveIncrement::AdaptiveIncrement(Json::Value const& increment_data)
{
    this->parse_input(increment_data);
}

void AdaptiveIncrement::update_convergence_state(bool is_converged)
{
    double constexpr cutback_factor = 0.5;
    double constexpr forward_factor = 2.0;

    if (is_converged)
    {
        consecutive_unconverged = 0;
        successful_increments++;

        double const factor_increment = std::min(forward_factor * current_factor, maximum_increment);

        is_applied = std::abs(current_factor - final_factor) < 1.0e-10 &&
                     current_factor + factor_increment > final_factor;

        if (is_applied) return;

        known_good_factor = current_factor;

        current_factor = std::min(factor_increment + current_factor, final_factor);

        std::cout << termcolor::green
                  << "\n****************************************************************\n\n"
                  << "  Convergence detected - setting factor to " << current_factor
                  << " for next attempt\n"
                  << "\n****************************************************************\n\n"
                  << termcolor::reset;
        return;
    }

    double const factor_increment =
        std::max(minimum_increment, (current_factor - known_good_factor) * cutback_factor);

    current_factor = known_good_factor + factor_increment;

    std::cout << termcolor::yellow
              << "\n******************************************************************\n\n"
              << "  Non-convergence detected - adaptive cut back to " << current_factor
              << " performed.\n"
              << "\n******************************************************************\n\n"
              << termcolor::reset;

    consecutive_unconverged++;

    if (consecutive_unconverged == increment_limit)
    {
        throw std::runtime_error("Convergence not judged likely\n");
    }
}

void AdaptiveIncrement::reset(Json::Value const& new_increment_data)
{
    this->parse_input(new_increment_data);

    is_applied = false;
    known_good_factor = 0.0;
    consecutive_unconverged = 0;
}

void AdaptiveIncrement::parse_input(Json::Value const& increment_data)
{
    this->check_increment_data(increment_data);

    // Initial factor determined by input
    initial_factor = increment_data["Increment"]["Initial"].asDouble();

    minimum_increment = increment_data["Increment"]["Minimum"].asDouble();
    maximum_increment = increment_data["Increment"]["Maximum"].asDouble();

    is_ramp_variation = increment_data["LoadVariation"].asString() == "Ramp";

    // Final factor should be load fully applied
    final_factor = increment_data["Period"].asDouble();

    current_factor = initial_factor;
}

void AdaptiveIncrement::check_increment_data(Json::Value const& increment_data)
{
    if (increment_data["Period"].empty())
        throw std::runtime_error("Time data requires a \"Period\" value\n");

    if (increment_data["Increment"].empty())
        throw std::runtime_error("Increment data not provided!\n");

    if (increment_data["Increment"]["Initial"].empty())
        throw std::runtime_error("Increment-Initial data not provided!\n");

    if (increment_data["Increment"]["Minimum"].empty())
        throw std::runtime_error("Increment-Minimum data not provided!\n");

    if (increment_data["Increment"]["Maximum"].empty())
        throw std::runtime_error("Increment-Maximum data not provided!\n");

    if (increment_data["LoadVariation"].empty())
        throw std::runtime_error(
            "\"LoadVariation\" is not specified as \"Ramp\" or \"Instantaneous\"\n");
}
}
