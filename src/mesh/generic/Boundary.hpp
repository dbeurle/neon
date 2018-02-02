
#pragma once

#include <utility>
#include <vector>

#include "io/json_forward.hpp"

namespace neon
{
/**
 * Boundary is a base class for boundary conditions which performs the load
 * interpolation logic.
 */
class Boundary
{
public:
    explicit Boundary(json const& times, json const& loads);

    ~Boundary() = default;

    [[nodiscard]] std::vector<double> time_history() const;

    /**
     * Interpolates linearly between the old prescribed boundary value and
     * the newly prescribed boundary value.  This takes into account if the
     * load is ramped or propogated from the previous load/time step
     * @return an interpolated value
     */
    [[nodiscard]] double interpolate_prescribed_load(double const step_time) const;

protected:
    [[nodiscard]] double interpolate_prescribed_load(
        std::vector<std::pair<double, double>> const& time_value, double const step_time) const;

private:
    std::vector<std::pair<double, double>> time_load;
};
}
