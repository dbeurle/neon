
#pragma once

#include <tuple>
#include <vector>

#include <json/forwards.h>

namespace neon
{
class AdaptiveIncrement
{
public:
    AdaptiveIncrement(Json::Value const& increment_data);

    /** Check if the load increment is finalized */
    bool is_fully_applied() const { return is_applied; }

    /** Get the current scaling factor for the attempted increment */
    double factor() const { return current_factor; }

    /** Get the fraction through the adaptive incrementation */
    double load_factor() const { return is_ramp_variation ? current_factor / final_factor : 1.0; }

    /** Get the pseudo time step size */
    double increment() const { return current_factor - known_good_factor; }

    /** The number of steps taken so far in the incrementor */
    auto step() const { return successful_increments; }

    /** Update the convergence state to determine the next increment */
    void update_convergence_state(bool is_converged);

    void reset(Json::Value const& new_increment_data);

protected:
    void parse_input(Json::Value const& increment_data);

    void check_increment_data(Json::Value const& increment_data);

protected:
    int const increment_limit = 5; //!< Maximum allowable increments
    int successful_increments = 0;

    double initial_factor = 1.0;
    double final_factor = 1.0;
    double current_factor = 1.0;

    double minimum_increment;
    double maximum_increment;

    double known_good_factor = 0.0;

    bool is_applied = false;
    bool is_ramp_variation = true;

    int consecutive_unconverged = 0; //!< Count of consecutive unsuccessful attempts
};
}
