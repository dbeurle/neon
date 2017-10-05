
#include "Boundary.hpp"

#include "numeric/DenseTypes.hpp"

#include <functional>

#include <range/v3/algorithm.hpp>
#include <range/v3/view.hpp>

#include <json/value.h>

namespace neon
{
Boundary::Boundary(Json::Value const& times, Json::Value const& loads)
{
    using namespace ranges;

    if (times.size() < 2)
    {
        throw std::runtime_error("\"Time\" is not a vector of at least two elements");
    }

    if (loads.size() < 2)
    {
        throw std::runtime_error("Boundary value is not a vector of at least two elements");
    }

    if (times.size() != loads.size())
    {
        throw std::runtime_error("The time and load vectors must be of the same size.");
    }
    if (ranges::adjacent_find(times,
                              [](auto const& left_value, auto const& right_value) {
                                  return left_value.asDouble() > right_value.asDouble();
                              })
        != times.end())
    {
        throw std::runtime_error("Load time must be monotonically increasing.");
    }

    time_load = view::zip(times, loads) | view::transform([](auto const& time_load) {
                    return std::make_pair(time_load.first.asDouble(), time_load.second.asDouble());
                });
}

std::vector<double> Boundary::time_history() const { return time_load | ranges::view::keys; }

double Boundary::interpolate_prescribed_load(double const step_time) const
{
    using ranges::adjacent_find;
    using ranges::find_if;
    using ranges::next;

    // Find if we have a time that matches exactly to a load
    if (auto match = find_if(time_load,
                             [&](auto const& value) { return is_approx(value.first, step_time); });
        match != time_load.end())
    {
        return (*match).second;
    }

    // Otherwise we need to interpolate between the values
    auto const lower_position = adjacent_find(time_load, [&](auto const& left, auto const& right) {
        return left.first < step_time && step_time < right.first;
    });

    auto const lower_pair = *lower_position;
    auto const upper_pair = *(ranges::next(lower_position));

    auto const[time_0, load_0] = lower_pair;
    auto const[time_1, load_1] = upper_pair;

    auto const interpolation_factor = (step_time - time_0) / (time_1 - time_0);

    return (1.0 - interpolation_factor) * load_0 + interpolation_factor * load_1;
}
}
