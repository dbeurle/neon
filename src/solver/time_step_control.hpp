
#pragma once

#include "io/json_forward.hpp"
#include "numeric/index_types.hpp"

namespace neon
{
class time_step_control
{
public:
    /** Constructor for static solution */
    time_step_control() = default;

    /** Constructor for dynamic solution */
    time_step_control(json const& time_data);

    [[nodiscard]] double current_time_step_size() const { return time_step_size; }

    [[nodiscard]] std::int64_t number_of_time_steps() const { return time_steps; }

    void increment() { ++current_time_step; }

    bool is_finished() const { return current_time_step == time_steps; }

protected:
    std::int64_t time_steps{1}, current_time_step{0};
    double time_step_size{1.0};
};
}
