
#include "AdaptiveLoadStep.hpp"

#include <exception>
#include <range/v3/algorithm.hpp>
#include <termcolor/termcolor.hpp>

#include <json/value.h>

namespace neon
{
AdaptiveLoadStep::AdaptiveLoadStep(Json::Value const& increment_data)
{
    this->parse_input(increment_data);
}

void AdaptiveLoadStep::update_convergence_state(bool is_converged)
{
    double constexpr cutback_factor = 0.5;
    double constexpr forward_factor = 2.0;

    int constexpr terminal_indent = 4;

    if (is_converged)
    {
        consecutive_unconverged = 0;
        successful_increments++;

        double const dt = std::min(forward_factor * current_time, maximum_increment);

        is_applied = std::abs(current_time - final_time) < 1.0e-10 &&
                     current_time + dt > final_time;

        last_converged_time = current_time;

        if (is_applied) return;

        current_time = std::min(dt + current_time, final_time);

        std::cout << termcolor::green << termcolor::bold
                  << std::string(terminal_indent, ' ')
                  << "Convergence detected - step time set to " << current_time
                  << " for next attempt\n"
                  << termcolor::reset << std::flush;
    }
    else
    {
        double const dt = std::max(minimum_increment,
                                   (current_time - last_converged_time) * cutback_factor);

        current_time = last_converged_time + dt;

        std::cout << termcolor::yellow << termcolor::bold
                  << std::string(terminal_indent, ' ')
                  << "\nNon-convergence detected - step time set to " << current_time
                  << "\n"
                  << termcolor::reset << std::flush;

        consecutive_unconverged++;

        if (consecutive_unconverged == increment_limit)
        {
            throw std::runtime_error(
                "Too many reductions.  Convergence not judged likely\n");
        }
    }
}

void AdaptiveLoadStep::reset(Json::Value const& new_increment_data)
{
    // Update the history counters
    total_time += current_time;

    this->parse_input(new_increment_data);

    is_applied = false;
    last_converged_time = 0.0;
    consecutive_unconverged = 0;
}

void AdaptiveLoadStep::parse_input(Json::Value const& increment_data)
{
    this->check_increment_data(increment_data);

    // Initial factor determined by input
    initial_time = increment_data["Increments"]["Initial"].asDouble();

    minimum_increment = increment_data["Increments"]["Minimum"].asDouble();
    maximum_increment = increment_data["Increments"]["Maximum"].asDouble();

    final_time = increment_data["Period"].asDouble();

    current_time = initial_time;
}

void AdaptiveLoadStep::check_increment_data(Json::Value const& increment_data)
{
    if (increment_data["Period"].empty())
        throw std::runtime_error("Time data requires a \"Period\" value\n");

    if (increment_data["Increments"].empty())
        throw std::runtime_error("\"Increments\" not provided!\n");

    if (increment_data["Increments"]["Initial"].empty())
        throw std::runtime_error("Increment-Initial data not provided!\n");

    if (increment_data["Increments"]["Minimum"].empty())
        throw std::runtime_error("Increment-Minimum data not provided!\n");

    if (increment_data["Increments"]["Maximum"].empty())
        throw std::runtime_error("Increment-Maximum data not provided!\n");
}
}
