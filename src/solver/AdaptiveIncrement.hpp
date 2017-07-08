
#pragma once

#include <tuple>
#include <vector>

namespace neon
{
class AdaptiveIncrement
{
public:
    AdaptiveIncrement(double initial_factor, double final_factor, int attemptLimit);

    /** Check if the load increment is finalized */
    bool is_fully_applied() const { return is_applied; }

    /** Get the current scaling factor for the attempted increment */
    double factor() const { return current_factor; }

    auto increment() const { return successful_increments; }

    /** Update the convergence state to determine the next increment */
    void update_convergence_state(bool is_converged);

protected:
    int increment_attempts = 0; //!< Total number of attempts at a given factor
    int increment_limit;        //!< Maximum allowable increments
    int successful_increments = 0;

    double initial_factor;
    double final_factor;
    double current_factor;

    double last_known_good_factor = 0.0;
    bool is_applied = false;

    double cutback_factor = 0.5; //!< How much to cut back by when not converged
    double forward_factor = 2.0; //!< How much to increase when converged

    int consecutive_unconverged = 0; //!< Count of consecutive unsuccessful attempts
};
}
