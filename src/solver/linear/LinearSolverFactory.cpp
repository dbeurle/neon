
#include "LinearSolverFactory.hpp"

#ifdef ENABLE_CUDA
#include "ConjugateGradientGPU.hpp"
#endif

#include "MUMPS.hpp"
#include "PaStiX.hpp"

#include <exception>

#include "io/json.hpp"

namespace neon
{
std::unique_ptr<LinearSolver> make_linear_solver(json const& solver_data,
                                                 bool const is_symmetric)
{
    std::string const& solver_name = solver_data["Solver"].asString();

    if (solver_name == "PaStiX")
    {
        if (is_symmetric) return std::make_unique<PaStiXLDLT>();
        return std::make_unique<PaStiXLU>();
    }
    else if (solver_name == "MUMPS")
    {
        if (is_symmetric) return std::make_unique<MUMPSLLT>();
        return std::make_unique<MUMPSLU>();
    }
    else if (solver_name == "Direct")
    {
        if (is_symmetric) return std::make_unique<SparseLLT>();

        return std::make_unique<SparseLU>();
    }
    else if (solver_name == "Iterative")
    {
        if (solver_data.isMember("Tolerance") && solver_data.isMember("MaxIterations"))
        {
            if (is_symmetric)
            {
                return std::make_unique<ConjugateGradient>(solver_data["Tolerance"].asDouble(),
                                                           solver_data["MaxIterations"].asInt());
            }
            return std::make_unique<BiCGStab>(solver_data["Tolerance"].asDouble(),
                                              solver_data["MaxIterations"].asInt());
        }
        else if (solver_data.isMember("Tolerance"))
        {
            if (is_symmetric)
            {
                return std::make_unique<ConjugateGradient>(solver_data["Tolerance"].asDouble());
            }
            return std::make_unique<BiCGStab>(solver_data["Tolerance"].asDouble());
        }
        else if (solver_data.isMember("MaxIterations"))
        {
            if (is_symmetric)
            {
                return std::make_unique<ConjugateGradient>(solver_data["MaxIterations"].asInt());
            }
            return std::make_unique<BiCGStab>(solver_data["MaxIterations"].asInt());
        }
        else
        {
            if (is_symmetric) return std::make_unique<ConjugateGradient>();
            return std::make_unique<BiCGStab>();
        }
    }

    else if (solver_name == "IterativeGPU")
    {
#ifdef ENABLE_CUDA

        if (!is_symmetric)
        {
            throw std::runtime_error("A non-symmetric iterative solver for GPU has not yet been "
                                     "implemented\n");
        }

        if (solver_data.isMember("Tolerance") && solver_data.isMember("MaxIterations"))
        {
            return std::make_unique<ConjugateGradientGPU>(solver_data["Tolerance"].asDouble(),
                                                          solver_data["MaxIterations"].asInt());
        }
        else if (solver_data.isMember("Tolerance"))
        {
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
                                 "configured with -DENABLE_CUDA=ON\n");
#endif
    }
    else
    {
        throw std::runtime_error("Did not find a linear solver\n");
    }
    return nullptr;
}
}
