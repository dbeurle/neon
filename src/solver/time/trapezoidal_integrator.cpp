
#include "trapezoidal_integrator.hpp"

#include <iostream>

#include "io/json.hpp"

namespace neon
{
trapezoidal_integrator::trapezoidal_integrator(json const& time_data)
{
    for (auto const f : {"Increments", "Period", "Method"})
    {
        if (!time_data.count(f))
        {
            throw std::runtime_error("The \"Time\" field requires a \"" + std::string(f) + "\" field");
        }
    }

    if (!time_data["Increments"].count("Initial"))
    {
        throw std::runtime_error("An \"Initial\" time field is required");
    }

    std::cout << "\n";

    if (auto const& solver = time_data["Method"].get<std::string>(); solver == "ExplicitEuler")
    {
        std::cout << std::string(4, ' ') << "Time discretisation is the explicit Euler method\n";
        method = 0.0;
    }
    else if (solver == "CrankNicolson")
    {
        std::cout << std::string(4, ' ') << "Time discretisation is the Crank-Nicolson method\n";
        method = 0.5;
    }
    else if (solver == "ImplicitEuler")
    {
        std::cout << std::string(4, ' ') << "Time discretisation is the implicit Euler method\n";
        method = 1.0;
    }
    else
    {
        throw std::runtime_error("Valid types of time solvers are \"ExplicitEuler\", "
                                 "\"ImplicitEuler\" or \"CrankNicolson\"");
    }
    // Read in the time data
    start_time = 0.0;
    final_time = time_data["Period"];
    time_step_size = time_data["Increments"]["Initial"];

    std::cout << std::string(4, ' ') << "Start time     : " << start_time << std::endl;
    std::cout << std::string(4, ' ') << "Final time     : " << final_time << std::endl;
    std::cout << std::string(4, ' ') << "Time step size : " << time_step_size << std::endl;
}

bool trapezoidal_integrator::loop()
{
    current_time_step++;

    time += time_step_size;

    return time < final_time;
}
}
