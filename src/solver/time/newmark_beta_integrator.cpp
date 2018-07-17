
#include "newmark_beta_integrator.hpp"

#include "io/json.hpp"
#include <iostream>

namespace neon
{
newmark_beta_integrator::newmark_beta_integrator(json const& time_solver_data)
    : time_control(time_solver_data)
{
    if (time_solver_data["IntegrationOptions"].empty())
    {
        // Choose sensible defaults (if they exist)
        artifical_viscosity = 0.5;
        beta_parameter = 0.25;
        std::cout << "Warning: Newmark-Beta parameters were not specified.\nThe "
                     "artifical viscosity is set to 0.5 and the beta parameter is set to "
                     "0.25.\nYou can suppress this warning by using \"ViscousDamping\" and "
                     "\"BetaParameter\" under the \"IntegrationOptions\" field.";
    }
    else
    {
        if (time_solver_data["IntegrationOptions"]["ViscousDamping"].is_null())
        {
            throw std::runtime_error("IterationOptions - ViscousDamping was not set\n");
        }
        if (time_solver_data["IntegrationOptions"]["BetaParameter"].is_null())
        {
            throw std::runtime_error("IterationOptions - BetaParameter was not set\n");
        }

        artifical_viscosity = time_solver_data["IntegrationOptions"]["ViscousDamping"];

        beta_parameter = time_solver_data["IntegrationOptions"]["BetaParameter"];
    }

    if (are_parameters_unstable())
    {
        throw std::runtime_error("Chosen Newmark-Beta parameters are not stable\n");
    }
}

bool newmark_beta_integrator::time_loop()
{
    time_control.increment();
    return !time_control.is_finished();
}

bool newmark_beta_integrator::are_parameters_unstable() const
{
    // Unconditional stability condition
    return beta_parameter < artifical_viscosity / 2.0 || artifical_viscosity / 2.0 < 0.25;
}

double newmark_beta_integrator::mass_scaling_factor() const
{
    auto const time_step_size = time_control.current_time_step_size();
    return 1.0 / (beta_parameter * std::pow(time_step_size, 2));
}

double newmark_beta_integrator::time_step_size() const
{
    return time_control.current_time_step_size();
}

vector newmark_beta_integrator::approximate_displacements(vector const& a,
                                                          vector const& v,
                                                          vector const& d) const
{
    auto const time_step_size = time_control.current_time_step_size();

    return d + time_step_size * v
           + std::pow(time_step_size, 2) / 2.0 * (1.0 - 2.0 * beta_parameter) * a;
}

vector newmark_beta_integrator::approximate_velocities(vector const& a, vector const& v) const
{
    auto const time_step_size = time_control.current_time_step_size();
    return v + (1.0 - artifical_viscosity) * time_step_size * a;
}

vector newmark_beta_integrator::accelerations(vector const& a, vector const& v, vector const& d)
{
    auto const time_step_size = time_control.current_time_step_size();
    a_old = a;

    return 1.0 / (beta_parameter * std::pow(time_step_size, 2))
           * (d - approximate_displacements(a, v, d));
}
vector newmark_beta_integrator::velocities(vector const& a, vector const& v) const
{
    auto const time_step_size = time_control.current_time_step_size();
    return approximate_velocities(a_old, v) + artifical_viscosity * time_step_size * a;
}
}
