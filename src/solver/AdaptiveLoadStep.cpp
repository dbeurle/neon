
#include "AdaptiveLoadStep.hpp"

#include "numeric/DenseTypes.hpp"

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
    auto constexpr cutback_factor = 0.5;
    auto constexpr forward_factor = 2.0;

    auto constexpr terminal_indent = 4;

    if (is_converged)
    {
        consecutive_unconverged = 0;
        successful_increments++;
        last_converged_time = current_time;

        is_applied = is_approx(current_time, final_time);

        if (is_applied) return;

        double const new_time = std::min(forward_factor * current_time,
                                         current_time + maximum_increment);

        current_time = std::min(new_time, final_time);

        std::cout << std::string(terminal_indent, ' ') << termcolor::green
                  << termcolor::bold << "Convergence detected - step time set to "
                  << current_time << " for next attempt\n"
                  << termcolor::reset << std::flush;
    }
    else
    {
        double const dt = std::max(minimum_increment,
                                   (current_time - last_converged_time) * cutback_factor);

        current_time = last_converged_time + dt;

        if (current_time > final_time)
        {
            throw std::runtime_error("Minimum increment is not small enough to resolve "
                                     "the step\n");
        }

        if (current_time < total_time)
            throw std::runtime_error("Step has reduced below the final time for the "
                                     "last load case\n");

        std::cout << "\n"
                  << std::string(terminal_indent, ' ') << termcolor::yellow
                  << termcolor::bold << "Non-convergence detected - step time set to "
                  << current_time << "\n"
                  << termcolor::reset << std::flush;

        consecutive_unconverged++;

        if (consecutive_unconverged == increment_limit)
        {
            throw std::runtime_error("Too many reductions.  Convergence not judged "
                                     "likely\n");
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
