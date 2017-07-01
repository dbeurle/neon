
#pragma once

#include "LinearSolver.hpp"

#include <iostream>
#include <json/json.h>
#include <memory>

namespace neon
{
class LinearSolverFactory
{
public:
    static std::unique_ptr<LinearSolver> make(Json::Value const& solver_data)
    {
        std::string const& solver_name = solver_data["Solver"].asString();
        if (solver_name == "Pastix")
        {
            return std::make_unique<PaStiX>();
        }
        else if (solver_name == "MUMPS")
        {
            return std::make_unique<MUMPS>();
        }
        else if (solver_name == "SparseLU")
        {
            return std::make_unique<SparseLU>();
        }
        else if (solver_name == "pCG")
        {
            if (not solver_data["Tolerance"].empty() and not solver_data["MaxIterations"].empty())
            {
                return std::make_unique<pCG>(solver_data["Tolerance"].asDouble(),
                                             solver_data["MaxIterations"].asInt());
            }
            else if (not solver_data["Tolerance"].empty())
            {
                // todo add messages in these statements to print out
                // to a file that other options have not been set
                // and explain how to set them
                return std::make_unique<pCG>(solver_data["Tolerance"].asDouble());
            }
            else if (not solver_data["MaxIterations"].empty())
            {
                return std::make_unique<pCG>(solver_data["MaxIterations"].asInt());
            }
            else
            {
                return std::make_unique<pCG>();
            }
        }
        else if (solver_name == "BiCGStab")
        {
            if (not solver_data["Tolerance"].empty() and not solver_data["MaxIterations"].empty())
            {
                return std::make_unique<BiCGSTAB>(solver_data["Tolerance"].asDouble(),
                                                  solver_data["MaxIterations"].asInt());
            }
            else if (not solver_data["Tolerance"].empty())
            {
                // todo add messages in these statements to print out
                // to a file that other options have not been set
                // and explain how to set them
                return std::make_unique<BiCGSTAB>(solver_data["Tolerance"].asDouble());
            }
            else if (!solver_data["MaxIterations"].empty())
            {
                return std::make_unique<BiCGSTAB>(solver_data["MaxIterations"].asInt());
            }
            else
            {
                return std::make_unique<BiCGSTAB>();
            }
        }
        else
        {
            std::cerr << "Did not find a linear solver\n";
            // TODO exceptions ;)
        }
        return nullptr;
    }
};
}
