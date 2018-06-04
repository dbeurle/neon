
#pragma once

#include "io/json_forward.hpp"

#include <utility>
#include <vector>

namespace neon
{
/**
 * boundary is a base class for boundary conditions which performs the load
 * interpolation logic.
 */
class boundary
{
public:
    explicit boundary(json const& times, json const& loads);

    explicit boundary(json const& boundary, std::string const& name, double const generate_time_step);

    ~boundary() = default;

    /// \return true if the boundary is active \sa is_not_active
    [[nodiscard]] bool is_active(double const time) const;

    /// \return true if the boundary is inactive \sa is_active
    [[nodiscard]] bool is_not_active(double const time) const;

    [[nodiscard]] std::vector<double> time_history() const;

    /// Interpolates linearly between the old prescribed boundary value and
    /// the newly prescribed boundary value.  This takes into account if the
    /// load is ramped or propogated from the previous load/time step
    /// \return an interpolated value
    [[nodiscard]] double interpolate_prescribed_load(double const step_time) const;

protected:
    [[nodiscard]] double interpolate_prescribed_load(
        std::vector<std::pair<double, double>> const& time_value,
        double const step_time) const;

    /// Allocate and check the time/load boundary values
    void allocate_time_load(json const& times, json const& loads);

    /**
     * Generate boundary values that follows a sinusoidal wave with fixed time step.
     * @param JSON boundary object
     * @param dofs {x,y,z}
     * @param time step
     */
    void generate_sinusoidal(json const& boundary,
                             std::string const& name,
                             double const generate_time_step);

private:
    std::vector<std::pair<double, double>> time_load;
};
}
