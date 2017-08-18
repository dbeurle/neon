
#include "LinearSolver.hpp"

#ifdef ENABLE_CUDA
#include "ConjugateGradientGPU.hpp"
#endif

#include "MUMPS.hpp"
#include "Pastix.hpp"

#include <exception>
#include <memory>
#include <string>

#include <json/json.h>

namespace neon
{
std::unique_ptr<LinearSolver> make_linear_solver(Json::Value const& solver_data)
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
    else if (solver_name == "ConjugateGradient")
    {
        if (solver_data.isMember("Tolerance") && solver_data.isMember("MaxIterations"))
        {
            return std::make_unique<ConjugateGradient>(solver_data["Tolerance"].asDouble(),
                                                       solver_data["MaxIterations"].asInt());
        }
        else if (solver_data.isMember("Tolerance"))
        {
            // todo add messages in these statements to print out
            // to a file that other options have not been set
            // and explain how to set them
            return std::make_unique<ConjugateGradient>(solver_data["Tolerance"].asDouble());
        }
        else if (solver_data.isMember("MaxIterations"))
        {
            return std::make_unique<ConjugateGradient>(solver_data["MaxIterations"].asInt());
        }
        else
        {
            return std::make_unique<ConjugateGradient>();
        }
    }

    else if (solver_name == "ConjugateGradientGPU")
    {
#ifdef ENABLE_CUDA
        if (solver_data.isMember("Tolerance") && solver_data.isMember("MaxIterations"))
        {
            return std::make_unique<ConjugateGradientGPU>(solver_data["Tolerance"].asDouble(),
                                                          solver_data["MaxIterations"].asInt());
        }
        else if (solver_data.isMember("Tolerance"))
        {
            // todo add messages in these statements to print out
            // to a file that other options have not been set
            // and explain how to set them
            return std::make_unique<ConjugateGradientGPU>(solver_data["Tolerance"].asDouble());
        }
        else if (solver_data.isMember("MaxIterations"))
        {
            return std::make_unique<ConjugateGradientGPU>(solver_data["MaxIterations"].asInt());
        }
        else
        {
            return std::make_unique<ConjugateGradientGPU>();
        }
#else
        throw std::runtime_error("ConjugateGradientGPU is only available when neon is "
                                 "compiled with -DENABLE_CUDA=ON\n");
#endif
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
        throw std::runtime_error("Did not find a linear solver\n");
    }
    return nullptr;
}
}
