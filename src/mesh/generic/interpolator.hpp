
#pragma once

#include <utility>
#include <vector>

#include <json/forwards.h>

namespace neon::boundary
{
/**
 * base is a base class for boundary conditions which performs the load
 * interpolation logic.
 */
class interpolator
{
public:
    explicit interpolator(Json::Value const& times, Json::Value const& loads);

    ~interpolator() = default;

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
