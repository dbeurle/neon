
#include "AdaptiveIncrement.hpp"

#include <exception>
#include <range/v3/algorithm.hpp>
#include <termcolor/termcolor.hpp>

namespace neon
{
AdaptiveIncrement::AdaptiveIncrement(double initial_factor, double final_factor, int increment_limit)
    : increment_limit(increment_limit),
      initial_factor(initial_factor),
      final_factor(final_factor),
      current_factor(initial_factor)
{
}

void AdaptiveIncrement::update_convergence_state(bool is_converged)
{
    if (is_converged)
    {
        is_applied = current_factor * forward_factor > final_factor &&
                     std::abs(final_factor - current_factor) < 1.0e-10;

        if (is_applied) return;

        last_known_good_factor = current_factor;

        current_factor = std::min(current_factor * forward_factor, final_factor);

        consecutive_unconverged = 0;

        successful_increments++;

        std::cout << termcolor::green
                  << "\n****************************************************************\n\n"
                  << "  Convergence detected - setting factor to " << current_factor
                  << " for next attempt\n"
                  << "\n****************************************************************\n\n"
                  << termcolor::reset;
    }
    else
    {
        current_factor = last_known_good_factor + (current_factor - last_known_good_factor) / 2.0;

        std::cout << termcolor::yellow
                  << "\n******************************************************************\n\n"
                  << "  Non-convergence detected - adaptive cut back to " << current_factor
                  << " performed.\n"
                  << "\n******************************************************************\n\n"
                  << termcolor::reset;

        consecutive_unconverged++;

        if (consecutive_unconverged == 5)
        {
            throw std::runtime_error("Convergence not judged likely\n");
        }
    }

    if (increment_attempts >= increment_limit)
    {
        throw std::runtime_error("Maximum number of increments reached.  Increase the limit as "
                                 "convergence is likely.\n");
    }
    increment_attempts++;
}
}
