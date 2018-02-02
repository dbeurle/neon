
#pragma once

#include "io/json_forward.hpp"

namespace neon
{
class TimeControl
{
public:
    /** Constructor for static solution */
    TimeControl() = default;

    /** Constructor for dynamic solution */
    TimeControl(json const& time_data);

    [[nodiscard]] double current_time_step_size() const { return time_step_size; }

        [[nodiscard]] int32_t number_of_time_steps() const
    {
        return time_steps;
    }

    void increment() { ++current_time_step; }

    bool is_finished() const { return current_time_step == time_steps; }

protected:
    int32_t time_steps{1}, current_time_step{0};
    double time_step_size = 1.0;
};
}
