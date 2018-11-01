
#include "time_step_control.hpp"

#include "io/json.hpp"

#include <stdexcept>

namespace neon
{
time_step_control::time_step_control(json const& time_data)
{
    if (time_data.find("start") == time_data.end())
    {
        throw std::domain_error("start time not specified in input (\"start\")\n");
    }
    if (time_data.find("end") == time_data.end())
    {
        throw std::domain_error("end time not specified in input (\"end\")\n");
    }
    if (time_data.find("step_size") == time_data.end())
    {
        throw std::domain_error("\"step_size\" not specified in input\n");
    }

    auto const total_time = time_data["end"].get<double>() - time_data["start"].get<double>();

    time_step_size = time_data["step_size"];

    time_steps = static_cast<std::int64_t>(total_time / time_step_size);
}
}
