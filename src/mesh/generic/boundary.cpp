
#include "boundary.hpp"

#include "numeric/float_compare.hpp"

#include "io/json.hpp"

#include <range/v3/algorithm/adjacent_find.hpp>
#include <range/v3/algorithm/find_if.hpp>
#include <range/v3/view/map.hpp>
#include <range/v3/view/transform.hpp>
#include <range/v3/view/zip.hpp>
#include <range/v3/view/map.hpp>

#include <cmath>
#include <numeric>

namespace neon
{
boundary::boundary(json const& times, json const& loads) { allocate_time_load(times, loads); }

boundary::boundary(json const& boundary, std::string const& name, double const generate_time_step)
{
    if (boundary.find("GenerateType") != boundary.end())
    {
        if (auto const& type = boundary["GenerateType"].get<std::string>(); type == "Sinusoidal")
        {
            generate_sinusoidal(boundary, name, generate_time_step);
        }
    }
    else
    {
        allocate_time_load(boundary["Time"], boundary[name]);
    }
}

std::vector<double> boundary::time_history() const { return time_load | ranges::view::keys; }

double boundary::interpolate_prescribed_load(double const step_time) const
{
    return interpolate_prescribed_load(time_load, step_time);
}

double boundary::interpolate_prescribed_load(std::vector<std::pair<double, double>> const& time_value,
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

void boundary::generate_sinusoidal(json const& boundary,
                                   std::string const& name,
                                   double const generate_time_step)
{
    using namespace ranges;

    for (auto const& mandatory_field : {"Period", "Phase", "NumberOfCycles"})
    {
        if (!boundary.count(mandatory_field))
        {
            throw std::domain_error("\"" + std::string(mandatory_field)
                                    + "\" was not specified in \"BoundaryCondition\".");
        }
    }

    std::vector<double> const time = boundary["Time"];
    std::vector<double> const amplitude = boundary[name];
    std::vector<double> const period = boundary["Period"];
    std::vector<double> const number_of_cycles = boundary["NumberOfCycles"];
    std::vector<double> const phase = boundary["Phase"];

    if (!((amplitude.size() == period.size()) && (period.size() == phase.size())
          && (phase.size() == number_of_cycles.size())))
    {
        throw std::domain_error("The amplitude, period, number_of_cycles, time_step and phase "
                                "vectors must be of the same size.");
    }

    std::vector<double> times, loads;
    times.push_back(0);

    // Loop over the blocks
    for (std::size_t i = 0; i < amplitude.size(); ++i)
    {
        auto const omega = 2.0 * M_PI / period[i];

        auto const number_of_steps = number_of_cycles[i] * period[i] / generate_time_step;

        std::vector<double> time_block(number_of_steps);
        std::iota(time_block.begin(), time_block.end(), 1);

        auto const scale_factor = number_of_cycles[i] * period[i] / time_block.back();

        time_block = time_block | view::transform([&scale_factor, &times](double x) {
                         return x * scale_factor + times.back();
                     });

        times.insert(times.end(), time_block.begin(), time_block.end());

        loads = times | view::transform([i, &amplitude, omega, &phase](double x) {
                    return amplitude[i] * std::sin(omega * x + phase[i]);
                });
    }
    allocate_time_load(times, loads);
}

void boundary::allocate_time_load(json const& times, json const& loads)
{
    using namespace ranges;

    if (times.size() < 2)
    {
        throw std::domain_error("\"Time\" is not a vector of at least two elements");
    }

    if (loads.size() < 2)
    {
        throw std::domain_error("boundary value is not a vector of at least two elements");
    }

    if (times.size() != loads.size())
    {
        throw std::domain_error("The time and load vectors must be of the same size.");
    }
    if (!std::is_sorted(begin(times), end(times)))
    {
        throw std::domain_error("Load times must be monotonically increasing.");
    }

    time_load = view::zip(times, loads) | view::transform([](auto const time_load) {
                    return std::make_pair(time_load.first, time_load.second);
                });
}
}
