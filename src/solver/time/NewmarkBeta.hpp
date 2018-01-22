
#pragma once

#include "numeric/DenseMatrix.hpp"
#include "solver/TimeControl.hpp"

#include "io/json.hpp"

namespace neon
{
class NewmarkBeta
{
public:
    NewmarkBeta(json const& time_solver_data);

    [[nodiscard]] double mass_scaling_factor() const;

    [[nodiscard]] vector approximate_displacements(vector const& a,
                                                   vector const& v,
                                                   vector const& d) const;

    [[nodiscard]] vector approximate_velocities(vector const& a, vector const& v) const;

    [[nodiscard]] vector accelerations(vector const& a, vector const& v, vector const& d);

    [[nodiscard]] vector velocities(vector const& a, vector const& v) const;

    /** Perform the time integration until returns false */
    [[nodiscard]] bool time_loop();

    [[nodiscard]] double time_step_size() const { return time_control.current_time_step_size(); }

    protected : [[nodiscard]] bool are_parameters_unstable() const;

protected:
    TimeControl time_control;
    double artifical_viscosity;
    double beta_parameter;

    vector a_old; //!< Old acceleration
};

inline double NewmarkBeta::mass_scaling_factor() const
{
    auto const time_step_size = time_control.current_time_step_size();
    return 1.0 / (beta_parameter * std::pow(time_step_size, 2));
}

inline vector NewmarkBeta::approximate_displacements(vector const& a,
                                                     vector const& v,
                                                     vector const& d) const
{
    auto const time_step_size = time_control.current_time_step_size();

    return d + time_step_size * v
           + std::pow(time_step_size, 2) / 2.0 * (1.0 - 2.0 * beta_parameter) * a;
}

inline vector NewmarkBeta::approximate_velocities(vector const& a, vector const& v) const
{
    auto const time_step_size = time_control.current_time_step_size();
    return v + (1.0 - artifical_viscosity) * time_step_size * a;
}

inline vector NewmarkBeta::accelerations(vector const& a, vector const& v, vector const& d)
{
    auto const time_step_size = time_control.current_time_step_size();
    a_old = a;

    return 1.0 / (beta_parameter * std::pow(time_step_size, 2))
           * (d - approximate_displacements(a, v, d));
}
inline vector NewmarkBeta::velocities(vector const& a, vector const& v) const
{
    auto const time_step_size = time_control.current_time_step_size();
    return approximate_velocities(a_old, v) + artifical_viscosity * time_step_size * a;
}
}
