
#pragma once

#include <utility>
#include <vector>

#include <json/forwards.h>

namespace neon
{
/**
 * Boundary is a base class for boundary conditions which performs the
 * interpolation logic.  This class offers restart functionality when performing
 * multiple time steps or load steps in a non-linear simulation.

 * \sa Dirichlet
 * \sa NonFollowerLoad
 */
class Boundary
{
public:
    explicit Boundary(Json::Value const& times, Json::Value const& loads);

    /** Update the prescribed value and load application to ramped or instantaneous */
    void internal_restart(double const prescribed_value_new, bool const is_load_ramped = true);

    /** Maintains the previous load and does not ramp */
    void internal_restart();

    std::vector<double> time_history() const;

    /**
     * Interpolates linearly between the old prescribed boundary value and
     * the newly prescribed boundary value.  This takes into account if the
     * load is ramped or propogated from the previous load/time step
     * @return an interpolated value
     */
    double interpolate_prescribed_load(double const step_time) const;

private:
    std::vector<std::pair<double, double>> time_load;
};
}
