
#include "linear_solver_factory.hpp"

#ifdef ENABLE_CUDA
#include "conjugate_gradientGPU.hpp"
#endif

#include "MUMPS.hpp"
#include "PaStiX.hpp"
#include "io/json.hpp"

#include <exception>

namespace neon
{
std::unique_ptr<LinearSolver> make_linear_solver(json const& solver_data, bool const is_symmetric)
{
    std::string const& solver_name = solver_data["Type"].get<std::string>();

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
        if (solver_data.count("Tolerance") && solver_data.count("MaxIterations"))
        {
            if (is_symmetric)
            {
                return std::make_unique<conjugate_gradient>(solver_data["Tolerance"].get<double>(),
                                                            solver_data["MaxIterations"].get<int>());
            }
            return std::make_unique<
                biconjugate_gradient_stabilised>(solver_data["Tolerance"].get<double>(),
                                                 solver_data["MaxIterations"].get<int>());
        }
        else if (solver_data.count("Tolerance"))
        {
            if (is_symmetric)
            {
                return std::make_unique<conjugate_gradient>(solver_data["Tolerance"].get<double>());
            }
            return std::make_unique<biconjugate_gradient_stabilised>(
                solver_data["Tolerance"].get<double>());
        }
        else if (solver_data.count("MaxIterations"))
        {
            if (is_symmetric)
            {
                return std::make_unique<conjugate_gradient>(solver_data["MaxIterations"].get<int>());
            }
            return std::make_unique<biconjugate_gradient_stabilised>(
                solver_data["MaxIterations"].get<int>());
        }
        else
        {
            if (is_symmetric) return std::make_unique<conjugate_gradient>();
            return std::make_unique<biconjugate_gradient_stabilised>();
        }
    }

    else if (solver_name == "IterativeGPU")
    {
#ifdef ENABLE_CUDA

        if (!is_symmetric)
        {
            throw std::domain_error("A non-symmetric iterative solver for GPU has not yet been "
                                    "implemented\n");
        }

        if (solver_data.count("Tolerance") && solver_data.count("MaxIterations"))
        {
            return std::make_unique<conjugate_gradientGPU>(solver_data["Tolerance"].get<double>(),
                                                           solver_data["MaxIterations"].get<int>());
        }
        else if (solver_data.count("Tolerance"))
        {
            return std::make_unique<conjugate_gradientGPU>(solver_data["Tolerance"].get<double>());
        }
        else if (solver_data.count("MaxIterations"))
        {
            return std::make_unique<conjugate_gradientGPU>(solver_data["MaxIterations"].get<int>());
        }
        else
        {
            return std::make_unique<conjugate_gradientGPU>();
        }
#else
        throw std::domain_error("conjugate_gradientGPU is only available when neon is "
                                "configured with -DENABLE_CUDA=ON\n");
#endif
    }
    else
    {
        throw std::domain_error("Did not find a linear solver type.  Did you try specifying "
                                "\"Type\" for the linear solver name?\n");
    }
    return nullptr;
}
}
