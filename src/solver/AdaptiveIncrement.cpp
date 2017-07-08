
#include "AdaptiveIncrement.hpp"

#include <iostream>

#include <range/v3/algorithm.hpp>

namespace neon
{
AdaptiveIncrement::AdaptiveIncrement(double initial_factor, double final_factor, int increment_limit)
    : increment_limit(increment_limit),
      initial_factor(initial_factor),
      final_factor(final_factor),
      current_factor(initial_factor)
{
}

bool AdaptiveIncrement::is_fully_applied() const { return is_applied; }

double AdaptiveIncrement::factor() const { return current_factor; }

void AdaptiveIncrement::update_convergence_state(bool is_converged)
{
    if (is_converged)
    {
        is_applied = current_factor * forward_factor >= final_factor;

        last_known_good_factor = current_factor;

        current_factor = std::min(current_factor * forward_factor, final_factor);

        consecutive_unconverged = 0;
    }
    else
    {
        current_factor = last_known_good_factor + (current_factor - last_known_good_factor) / 2.0;

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
