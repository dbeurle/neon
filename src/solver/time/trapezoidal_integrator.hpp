
#pragma once

#include "io/json.hpp"

namespace neon
{
class trapezoidal_integrator
{
public:
    trapezoidal_integrator(json const& time_solver_data);

    /** Perform the time integration until returns false */
    bool loop();

    double current_time_step_size() const noexcept { return time_step_size; }

    int iteration() const noexcept { return current_time_step; }

    double current_time() const noexcept { return time; }

protected:
    double method; //!< 0.0 if forward Euler, 0.5 if Crank-Nicolson and 1.0 if backward Euler

    double start_time{0.0}, final_time{1.0}, time_step_size{1.0};

    double time{start_time};
    int current_time_step{0};
};
}
