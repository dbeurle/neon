
#include "trapezoidal_integrator.hpp"

#include <iostream>

#include "io/json.hpp"

namespace neon
{
trapezoidal_integrator::trapezoidal_integrator(json const& time_data)
{
    for (auto const f : {"increments", "period", "method"})
    {
        if (time_data.find(f) == end(time_data))
        {
            throw std::runtime_error("The \"time\" field requires a \"" + std::string(f) + "\" field");
        }
    }

    if (time_data["increments"].find("initial") == end(time_data["increments"]))
    {
        throw std::runtime_error("An \"initial\" time field is required");
    }

    std::cout << "\n";

    if (std::string const& solver = time_data["method"]; solver == "explicit_euler")
    {
        std::cout << std::string(4, ' ') << "Time discretisation is the explicit Euler method\n";
        method = 0.0;
    }
    else if (solver == "crank_nicolson")
    {
        std::cout << std::string(4, ' ') << "Time discretisation is the Crank-Nicolson method\n";
        method = 0.5;
    }
    else if (solver == "implicit_euler")
    {
        std::cout << std::string(4, ' ') << "Time discretisation is the implicit Euler method\n";
        method = 1.0;
    }
    else
    {
        throw std::runtime_error("Valid types of time solvers are \"explicit_euler\", "
                                 "\"implicit_euler\" or \"crank_nicolson\"");
    }
    // Read in the time data
    start_time = 0.0;
    final_time = time_data["period"];
    time_step_size = time_data["increments"]["initial"];

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
