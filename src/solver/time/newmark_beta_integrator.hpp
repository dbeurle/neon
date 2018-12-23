
#pragma once

/// @file

#include "numeric/dense_matrix.hpp"
#include "solver/time_step_control.hpp"

#include "io/json_forward.hpp"

namespace neon
{
class newmark_beta_integrator
{
public:
    newmark_beta_integrator(json const& time_solver_data);

    [[nodiscard]] double mass_scaling_factor() const;

    [[nodiscard]] vector approximate_displacements(vector const& a,
                                                   vector const& v,
                                                   vector const& d) const;

    [[nodiscard]] vector approximate_velocities(vector const& a, vector const& v) const;

    [[nodiscard]] vector accelerations(vector const& a, vector const& v, vector const& d);

    [[nodiscard]] vector velocities(vector const& a, vector const& v) const;

    /// Perform the time integration and return the state
    /// \return \p true if integration not finished, \p false otherwise.
    [[nodiscard]] bool time_loop();

    [[nodiscard]] double time_step_size() const;

protected:
    [[nodiscard]] bool are_parameters_unstable() const;

protected:
    time_step_control time_control;

    double artifical_viscosity;
    double beta_parameter;

    /// Old acceleration
    vector a_old;
};
}
