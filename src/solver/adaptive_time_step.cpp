
#include "adaptive_time_step.hpp"

#include "numeric/float_compare.hpp"
#include "io/json.hpp"

#include <exception>
#include <termcolor/termcolor.hpp>

namespace neon
{
adaptive_time_step::adaptive_time_step(json const& increment_data,
                                       std::vector<double> mandatory_time_history)
{
    // Set the time history to the input data
    for (auto const& t : mandatory_time_history) time_queue.push(t);

    // Remove the first entry if zero
    if (is_approx(time_queue.top(), 0.0)) time_queue.pop();

    parse_input(increment_data,
                *std::max_element(begin(mandatory_time_history), end(mandatory_time_history)));
}

void adaptive_time_step::update_convergence_state(bool const is_converged)
{
    auto constexpr cutback_factor{0.5};
    auto constexpr forward_factor{2.0};

    auto constexpr terminal_indent{4};

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
            throw std::domain_error("minimum increment is not small enough to resolve "
                                    "the time step\n");
        }
        if (current_time < 0.0)
        {
            throw std::domain_error("Step has reduced below the final time for the "
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
            throw std::domain_error("Too many reductions.  Convergence not judged "
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
    initial_time = increment_data["increments"]["initial"];

    bool const is_adaptive_increment = increment_data["increments"]["adaptive"];

    minimum_increment = is_adaptive_increment ? increment_data["increments"]["minimum"].get<double>()
                                              : initial_time;
    maximum_increment = is_adaptive_increment ? increment_data["increments"]["maximum"].get<double>()
                                              : initial_time;

    final_time = increment_data["period"];

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
    if (increment_data.find("period") == increment_data.end())
    {
        throw std::domain_error("Time data requires a \"period\" value\n");
    }
    if (increment_data.find("increments") == increment_data.end())
    {
        throw std::domain_error("\"increments\" not provided!\n");
    }

    auto const& increments_data = increment_data["increments"];

    if (increments_data.find("initial") == increments_data.end())
    {
        throw std::domain_error("Increment-initial data not provided!\n");
    }
    if (increments_data.find("minimum") == increments_data.end())
    {
        throw std::domain_error("Increment-minimum data not provided!\n");
    }
    if (increments_data.find("maximum") == increments_data.end())
    {
        throw std::domain_error("Increment-Maximum data not provided!\n");
    }
}

bool adaptive_time_step::is_highly_nonlinear() const
{
    return consecutive_unconverged > 0 || consecutive_converged < 4;
}
}
