
#include "Boundary.hpp"

#include "numeric/FloatingPointCompare.hpp"

#include <range/v3/algorithm/adjacent_find.hpp>
#include <range/v3/algorithm/find_if.hpp>

#include <range/v3/view/map.hpp>
#include <range/v3/view/transform.hpp>
#include <range/v3/view/zip.hpp>

#include "io/json.hpp"

namespace neon
{
Boundary::Boundary(json const& times, json const& loads)
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
        throw std::runtime_error("Load times must be monotonically increasing.");
    }

    time_load = view::zip(times, loads) | view::transform([](auto const& time_load) {
                    return std::make_pair(time_load.first.asDouble(), time_load.second.asDouble());
                });
}

std::vector<double> Boundary::time_history() const { return time_load | ranges::view::keys; }

double Boundary::interpolate_prescribed_load(double const step_time) const
{
    return interpolate_prescribed_load(time_load, step_time);
}

double Boundary::interpolate_prescribed_load(std::vector<std::pair<double, double>> const& time_value,
                                             double const step_time) const
{
    using ranges::adjacent_find;
    using ranges::find_if;
    using ranges::next;

    // Find if we have a time that matches exactly to a load
    if (auto match = find_if(time_value,
                             [&](auto const& value) { return is_approx(value.first, step_time); });
        match != time_value.end())
    {
        return (*match).second;
    }

    // Otherwise we need to interpolate between the values
    auto const lower_position = adjacent_find(time_value, [&](auto const& left, auto const& right) {
        return left.first < step_time && step_time < right.first;
    });

    auto const lower_pair = *lower_position;
    auto const upper_pair = *(ranges::next(lower_position));

    auto const [time_0, load_0] = lower_pair;
    auto const [time_1, load_1] = upper_pair;

    auto const interpolation_factor = (step_time - time_0) / (time_1 - time_0);

    return (1.0 - interpolation_factor) * load_0 + interpolation_factor * load_1;
}
}
