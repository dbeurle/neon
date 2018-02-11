
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
                     "artifical "
                     "viscosity is set to 0.5 and the beta parameter is set to "
                     "0.25.\nYou can "
                     "suppress this warning by using \"ViscousDamping\" and "
                     "\"BetaParameter\" "
                     "under the \"IntegrationOptions\" field.";
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
}
