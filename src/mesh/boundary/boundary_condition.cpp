
#include "boundary_condition.hpp"

#include "numeric/float_compare.hpp"
#include "math/linear_interpolation.hpp"

#include "io/json.hpp"

#include <range/v3/algorithm/adjacent_find.hpp>
#include <range/v3/algorithm/find_if.hpp>
#include <range/v3/view/map.hpp>
#include <range/v3/view/transform.hpp>
#include <range/v3/view/zip.hpp>

#include <cmath>
#include <numeric>

namespace neon
{
boundary_condition::boundary_condition(json const& times, json const& loads)
{
    allocate_time_load(times, loads);
}

boundary_condition::boundary_condition(json const& boundary_data,
                                       std::string const& name,
                                       double const generate_time_step)
{
    if (boundary_data.find("generate_type") != end(boundary_data))
    {
        if (std::string const& type = boundary_data["generate_type"]; type == "sinusoidal")
        {
            generate_sinusoidal(boundary_data, name, generate_time_step);
        }
    }
    else
    {
        allocate_time_load(boundary_data["time"], boundary_data[name]);
    }
}

bool boundary_condition::is_active(double const time) const
{
    return m_times.front() <= time && time <= m_times.back();
}

bool boundary_condition::is_not_active(double const time) const { return !is_active(time); }

auto boundary_condition::times() const -> std::vector<double> const& { return m_times; }

double boundary_condition::interpolate_prescribed_load(double const step_time) const
{
    return this->interpolate(m_times, m_values, step_time);
}

double boundary_condition::interpolate(std::vector<double> const& x_values,
                                       std::vector<double> const& y_values,
                                       double const x_value) const
{
    // Check if the time exact matches an existing time value
    if (auto match = std::find_if(begin(x_values),
                                  end(x_values),
                                  [&](auto const& value) { return is_approx(value, x_value); });
        match != end(x_values))
    {
        return y_values[std::distance(begin(x_values), match)];
    }

    // Otherwise we need to interpolate between the values
    auto [lower_value, upper_value] = std::equal_range(begin(x_values), end(x_values), x_value);

    // Handle edge case
    if (lower_value != begin(x_values))
    {
        lower_value = std::prev(lower_value);
    }

    if (lower_value == end(x_values) || upper_value == end(x_values) || lower_value == upper_value)
    {
        throw std::domain_error("step_time was not found as a valid time to interpolate");
    }
    return linear_interpolation(relative_distance(x_value, *lower_value, *upper_value),
                                y_values[std::distance(begin(x_values), lower_value)],
                                y_values[std::distance(begin(x_values), upper_value)]);
}

void boundary_condition::generate_sinusoidal(json const& boundary,
                                             std::string const& name,
                                             double const generate_time_step)
{
    using namespace ranges;

    for (auto const& mandatory_field : {"period", "phase", "number_of_cycles"})
    {
        if (boundary.find(mandatory_field) == end(boundary))
        {
            throw std::domain_error("\"" + std::string(mandatory_field)
                                    + "\" was not specified in \"boundaries\".");
        }
    }

    std::vector<double> const amplitude = boundary[name];
    std::vector<double> const period = boundary["period"];
    std::vector<double> const number_of_cycles = boundary["number_of_cycles"];
    std::vector<double> const phase = boundary["phase"];

    if (!(amplitude.size() == period.size() && period.size() == phase.size()
          && phase.size() == number_of_cycles.size()))
    {
        throw std::domain_error("The amplitude, period, number_of_cycles, time_step and phase "
                                "vectors must be all be the same size.");
    }

    // the first time step is zero
    std::vector<double> times(1, 0.0), loads(1, 0.0);

    // Loop over blocks of cycles.
    // Each block corresponds to one entry in amplitude, period, number_of_cycles and phase
    for (std::size_t i = 0; i < amplitude.size(); ++i)
    {
        auto const angular_frequency = 2.0 * M_PI / period[i];

        // number of time steps within one block of cycles
        auto const number_of_steps = number_of_cycles[i] * period[i] / generate_time_step;

        // steps over the current block excluding zero
        std::vector<double> time_block(number_of_steps);
        std::iota(begin(time_block), end(time_block), 1l);

        // scale the steps to time and shift if the current block is not the first
        auto const scale_factor = number_of_cycles[i] * period[i] / time_block.back();

        time_block = time_block | view::transform([&scale_factor, &times](double x) {
                         return x * scale_factor + times.back();
                     });

        // compute the load corresponding to each time step
        std::vector<double> load_block = time_block
                                         | view::transform([i, &amplitude, angular_frequency, &phase](
                                                               double x) {
                                               return amplitude[i]
                                                      * std::sin(angular_frequency * x + phase[i]);
                                           });

        // append the current block time steps and loads to the total time and load vectors
        times.insert(end(times), begin(time_block), end(time_block));
        loads.insert(end(loads), begin(load_block), end(load_block));
    }
    allocate_time_load(times, loads);
}

void boundary_condition::allocate_time_load(json const& times, json const& loads)
{
    m_times = times.get<std::vector<double>>();
    m_values = loads.get<std::vector<double>>();

    if (m_times.size() < 2)
    {
        throw std::domain_error("\"time\" is not a vector of at least two elements");
    }

    if (m_values.size() < 2)
    {
        throw std::domain_error("The specified boundary values do not contain at least two "
                                "elements");
    }

    if (m_times.size() != m_values.size())
    {
        throw std::domain_error("The time and load vectors must be of the same size.");
    }

    if (!std::is_sorted(begin(m_times), end(m_times)))
    {
        throw std::domain_error("Load times must be monotonically increasing.");
    }
}
}
