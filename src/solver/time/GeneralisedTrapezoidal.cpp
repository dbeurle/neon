
#include "GeneralisedTrapezoidal.hpp"

#include <iostream>

#include <json/value.h>

namespace neon
{
GeneralisedTrapezoidal::GeneralisedTrapezoidal(Json::Value const& time_data)
{
    for (auto const f : {"Increments", "Period", "Method"})
    {
        if (!time_data.isMember(f))
        {
            throw std::runtime_error("The \"Time\" field requires a \"" + std::string(f) + "\" field");
        }
    }

    if (!time_data["Increments"].isMember("Initial"))
    {
        throw std::runtime_error("An \"Initial\" time field is required");
    }

    if (auto const& solver = time_data["Method"].asString(); solver == "ExplicitEuler")
    {
        std::cout << "Using the explicit Euler method\n";
        method = 0.0;
    }
    else if (solver == "CrankNicolson")
    {
        std::cout << "Using the Crank-Nicolson method\n";
        method = 0.5;
    }
    else if (solver == "ImplicitEuler")
    {
        std::cout << "Using the implicit Euler method\n";
        method = 1.0;
    }
    else
    {
        throw std::runtime_error("Valid types of time solvers are \"ExplicitEuler\", "
                                 "\"ImplicitEuler\" or \"CrankNicolson\"");
    }
    // Read in the time data
    start_time = 0.0;
    final_time = time_data["Period"].asDouble();
    time_step_size = time_data["Increments"]["Initial"].asDouble();

    std::cout << "Start time = " << start_time << std::endl;
    std::cout << "Final time = " << final_time << std::endl;
    std::cout << "Time step size = " << time_step_size << std::endl;
}

bool GeneralisedTrapezoidal::loop()
{
    //
    return false;
}
}
