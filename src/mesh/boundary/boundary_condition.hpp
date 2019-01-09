
#pragma once

/// @file

#include "io/json_forward.hpp"

#include <utility>
#include <vector>

namespace neon
{
/// boundary is a base class for boundary conditions which performs the load
/// interpolation logic.
class boundary_condition
{
public:
    explicit boundary_condition(json const& times, json const& loads);

    explicit boundary_condition(json const& boundary_data,
                                std::string const& name,
                                double const generate_time_step);

    ~boundary_condition() = default;

    /// \return true if the boundary is active \sa is_not_active
    [[nodiscard]] bool is_active(double const time) const;

    /// \return true if the boundary is inactive \sa is_active
    [[nodiscard]] bool is_not_active(double const time) const;

    [[nodiscard]] auto times() const -> std::vector<double> const&;

    /// Interpolates linearly between the old prescribed boundary value and
    /// the newly prescribed boundary value.  This takes into account if the
    /// load is ramped or propogated from the previous load/time step
    /// \return an interpolated value
    [[nodiscard]] double interpolate_prescribed_load(double const step_time) const;

protected:
    [[nodiscard]] double interpolate(std::vector<double> const& x_values,
                                     std::vector<double> const& y_values,
                                     double const x_value) const;

    /// Allocate and check the time/load boundary values
    void allocate_time_load(json const& times, json const& loads);

    /// Generate boundary values that follows a sinusoidal wave with fixed time step. This may be
    /// used to generate block loading where each block follows a specific sinusoidal wave.
    /// \param boundary JSON boundary object
    /// \param Name of
    /// \param generate_time_step time step
    void generate_sinusoidal(json const& boundary,
                             std::string const& name,
                             double const generate_time_step);

protected:
    std::vector<double> m_times;
    std::vector<double> m_values;
};
}
