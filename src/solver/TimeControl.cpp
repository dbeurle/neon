
#include "TimeControl.hpp"

#include <exception>
#include "io/json.hpp"

namespace neon
{
TimeControl::TimeControl(json const& time_data)
{
    if (time_data["Start"].empty())
        throw std::runtime_error("Start time not specified in input (\"Start\")\n");

    if (time_data["End"].empty())
        throw std::runtime_error("End time not specified in input (\"End\")\n");

    if (time_data["StepSize"].empty())
        throw std::runtime_error("\"StepSize\" not specified in input\n");

    double total_time = time_data["End"].asDouble() - time_data["Start"].asDouble();

    time_step_size = time_data["StepSize"].asDouble();

    time_steps = static_cast<int>(total_time / time_step_size);
}
}
