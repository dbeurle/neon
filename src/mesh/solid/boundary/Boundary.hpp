
#pragma once

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
    explicit Boundary(double const prescribed_value, bool const is_load_ramped);

    /** Update the prescribed value and load application to ramped or instantaneous */
    void internal_restart(double const prescribed_value_new, bool const is_load_ramped = true);

    /** Maintains the previous load and does not ramp */
    void internal_restart();

    /**
     * Interpolates linearly between the old prescribed boundary value and
     * the newly prescribed boundary value.  This takes into account if the
     * load is ramped or propogated from the previous load/time step
     * @return an interpolated value
     */
    auto interpolate_prescribed_value(double const load_factor) const
    {
        return is_load_ramped ? (value_new - value_old) * load_factor + value_old : value_new;
    }

private:
    bool is_load_ramped;

    double value_old = 0.0;
    double value_new = 0.0;
};
}
