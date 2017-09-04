
#include "GeneralisedTrapezoidal.hpp"

#include <json/value.h>

namespace neon
{
GeneralisedTrapezoidal::GeneralisedTrapezoidal(Json::Value const& time_solver_data)
{
    if (auto const& solver = time_solver_data["Solver"].asString(); solver == "ExplicitEuler")
    {
        method = 0.0;
    }
    else if (solver == "CrankNicolson")
    {
        method = 0.5;
    }
    else if (solver == "BackwardEuler")
    {
        method = 1.0;
    }
    else
    {
        throw std::runtime_error("Valid types of hyperbolic solvers are \"ExplicitEuler\", "
                                 "\"BackwardEuler\" or \"CrankNicolson\"");
    }
}
}
