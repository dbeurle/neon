
#include "adaptive_time_step.hpp"

#include "numeric/float_compare.hpp"

#include <exception>
#include <termcolor/termcolor.hpp>

#include <range/v3/algorithm/max_element.hpp>

#include "io/json.hpp"

namespace neon
{
adaptive_time_step::adaptive_time_step(json const& increment_data,
                                   std::vector<double> mandatory_time_history)
{
    // Set the time history to the input data
    for (auto const& t : mandatory_time_history) time_queue.push(t);

    // Remove the first entry if it's zero
    if (is_approx(time_queue.top(), 0.0)) time_queue.pop();

    parse_input(increment_data, *ranges::max_element(mandatory_time_history));
}

void adaptive_time_step::update_convergence_state(bool const is_converged)
{
    auto constexpr cutback_factor = 0.5;
    auto constexpr forward_factor = 2.0;

    auto constexpr terminal_indent = 4;

    if (is_converged)
    {
        if (current_time >= time_queue.top()) time_queue.pop();

        last_converged_time_step_size = current_time - last_converged_time;
        last_converged_time = current_time;

        is_applied = is_approx(current_time, final_time) || time_queue.empty();

        consecutive_unconverged = 0;
        consecutive_converged++;
        successful_increments++;

        if (is_applied) return;

        // If the previous iterations required cut backs, then the next steps
        // should proceed slowly in the nonlinear region
        double const new_time = std::min(is_highly_nonlinear()
                                             ? last_converged_time_step_size + current_time
                                             : forward_factor * current_time,
                                         current_time + maximum_increment);

        current_time = std::min(time_queue.top(), std::min(new_time, final_time));

        std::cout << std::string(terminal_indent, ' ') << termcolor::green << termcolor::bold
                  << "Convergence detected - step time set to " << current_time
                  << " for next attempt\n"
                  << termcolor::reset << std::flush;
    }
    else
    {
        auto const dt = std::max(minimum_increment,
                                 (current_time - last_converged_time) * cutback_factor);

        current_time = last_converged_time + dt;

        if (current_time > final_time)
        {
            throw std::runtime_error("Minimum increment is not small enough to resolve "
                                     "the time step\n");
        }
        if (current_time < 0.0)
        {
            throw std::runtime_error("Step has reduced below the final time for the "
                                     "last load case\n");
        }

        std::cout << "\n"
                  << std::string(terminal_indent, ' ') << termcolor::yellow << termcolor::bold
                  << "Non-convergence detected - performing increment reduction.\n"
                  << termcolor::reset << std::flush;

        consecutive_unconverged++;
        consecutive_converged = 0;

        if (consecutive_unconverged == increment_limit)
        {
            throw std::runtime_error("Too many reductions.  Convergence not judged "
                                     "likely\n");
        }
    }
}

void adaptive_time_step::reset(json const& new_increment_data)
{
    // Update the history counters
    total_time += current_time;

    parse_input(new_increment_data, -1.0);

    is_applied = false;

    last_converged_time = last_converged_time_step_size = 0.0;
    consecutive_unconverged = consecutive_converged = 0;
}

void adaptive_time_step::parse_input(json const& increment_data, double const maximum_mandatory_time)
{
    check_increment_data(increment_data);

    // Initial factor determined by input
    initial_time = increment_data["Increments"]["Initial"];

    minimum_increment = increment_data["Increments"]["Minimum"];
    maximum_increment = increment_data["Increments"]["Maximum"];

    final_time = increment_data["Period"];

    if (maximum_mandatory_time > final_time)
    {
        std::cout << std::string(2, ' ') << termcolor::yellow << termcolor::bold
                  << "A boundary time is greater than the specified time period.  Please ensure "
                     "this is intended behaviour.\n"
                  << termcolor::reset << std::flush;
    }

    current_time = std::min(time_queue.top(), initial_time);
}

void adaptive_time_step::check_increment_data(json const& increment_data)
{
    if (!increment_data.count("Period"))
        throw std::runtime_error("Time data requires a \"Period\" value\n");

    if (!increment_data.count("Increments"))
        throw std::runtime_error("\"Increments\" not provided!\n");

    if (!increment_data["Increments"].count("Initial"))
        throw std::runtime_error("Increment-Initial data not provided!\n");

    if (!increment_data["Increments"].count("Minimum"))
        throw std::runtime_error("Increment-Minimum data not provided!\n");

    if (!increment_data["Increments"].count("Maximum"))
        throw std::runtime_error("Increment-Maximum data not provided!\n");
}
}
