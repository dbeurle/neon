
#include "time_step_control.hpp"

#include "io/json.hpp"

#include <stdexcept>

namespace neon
{
time_step_control::time_step_control(json const& time_data)
{
    if (time_data.find("Start") == time_data.end())
    {
        throw std::domain_error("Start time not specified in input (\"Start\")\n");
    }
    if (time_data.find("End") == time_data.end())
    {
        throw std::domain_error("End time not specified in input (\"End\")\n");
    }
    if (time_data.find("StepSize") == time_data.end())
    {
        throw std::domain_error("\"StepSize\" not specified in input\n");
    }

    auto const total_time = time_data["End"].get<double>() - time_data["Start"].get<double>();

    time_step_size = time_data["StepSize"];

    time_steps = static_cast<std::int64_t>(total_time / time_step_size);
}
}
