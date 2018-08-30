
#include "linear_solver_factory.hpp"

#ifdef ENABLE_CUDA
#include "conjugate_gradient_cuda.hpp"
#include "biconjugate_gradient_stabilised_cuda.hpp"
#endif

#ifdef ENABLE_OCL
#include "conjugate_gradient_ocl.hpp"
#include "biconjugate_gradient_stabilised_ocl.hpp"
#endif

#include "MUMPS.hpp"
#include "PaStiX.hpp"
#include "io/json.hpp"

#include <exception>

namespace neon
{
std::unique_ptr<linear_solver> make_linear_solver(json const& solver_data, bool const is_symmetric)
{
    std::string const& solver_name = solver_data["type"];

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
    else if (solver_name == "iterative")
    {
        if (solver_data.find("tolerance") != solver_data.end()
            && solver_data.find("maximum_iterations") != solver_data.end())
        {
            if (is_symmetric)
            {
                return std::make_unique<conjugate_gradient>(solver_data["tolerance"].get<double>(),
                                                            solver_data["maximum_iterations"]
                                                                .get<int>());
            }
            return std::make_unique<
                biconjugate_gradient_stabilised>(solver_data["tolerance"].get<double>(),
                                                 solver_data["maximum_iterations"].get<int>());
        }
        else if (solver_data.find("tolerance") != solver_data.end())
        {
            if (is_symmetric)
            {
                return std::make_unique<conjugate_gradient>(solver_data["tolerance"].get<double>());
            }
            return std::make_unique<biconjugate_gradient_stabilised>(
                solver_data["tolerance"].get<double>());
        }
        else if (solver_data.find("maximum_iterations") != solver_data.end())
        {
            if (is_symmetric)
            {
                return std::make_unique<conjugate_gradient>(
                    solver_data["maximum_iterations"].get<int>());
            }
            return std::make_unique<biconjugate_gradient_stabilised>(
                solver_data["maximum_iterations"].get<int>());
        }
        else
        {
            if (is_symmetric) return std::make_unique<conjugate_gradient>();
            return std::make_unique<biconjugate_gradient_stabilised>();
        }
    }

    else if (solver_name == "iterativeCUDA")
    {
#ifdef ENABLE_CUDA

        if (!is_symmetric)
        {
            if (solver_data.find("tolerance") != solver_data.end()
                && solver_data.find("maximum_iterations") != solver_data.end())
            {
                return std::make_unique<
                    biconjugate_gradient_stabilised_cuda>(solver_data["tolerance"].get<double>(),
                                                          solver_data["maximum_iterations"].get<int>());
            }
            else if (solver_data.find("tolerance") != solver_data.end())
            {
                return std::make_unique<biconjugate_gradient_stabilised_cuda>(
                    solver_data["tolerance"].get<double>());
            }
            else if (solver_data.find("maximum_iterations") != solver_data.end())
            {
                return std::make_unique<biconjugate_gradient_stabilised_cuda>(
                    solver_data["maximum_iterations"].get<int>());
            }
            return std::make_unique<biconjugate_gradient_stabilised_cuda>();
        }
        else
        {
            if (solver_data.find("tolerance") != solver_data.end()
                && solver_data.find("maximum_iterations") != solver_data.end())
            {
                return std::make_unique<conjugate_gradient_cuda>(solver_data["tolerance"].get<double>(),
                                                                 solver_data["maximum_iterations"]
                                                                     .get<int>());
            }
            else if (solver_data.find("tolerance") != solver_data.end())
            {
                return std::make_unique<conjugate_gradient_cuda>(
                    solver_data["tolerance"].get<double>());
            }
            else if (solver_data.find("maximum_iterations") != solver_data.end())
            {
                return std::make_unique<conjugate_gradient_cuda>(
                    solver_data["maximum_iterations"].get<int>());
            }
            return std::make_unique<conjugate_gradient_cuda>();
        }

#else
        throw std::domain_error("iterativeCUDA is only available when neon is "
                                "configured with -DENABLE_CUDA=1\n");
#endif
    }
    else if (solver_name == "iterativeOCL")
    {
#ifdef ENABLE_OCL
        if (!is_symmetric)
        {
            if (solver_data.find("tolerance") != solver_data.end()
                && solver_data.find("maximum_iterations") != solver_data.end())
            {
                return std::make_unique<
                    biconjugate_gradient_stabilised_ocl>(solver_data["tolerance"].get<double>(),
                                                         solver_data["maximum_iterations"].get<int>());
            }
            else if (solver_data.find("tolerance") != solver_data.end())
            {
                return std::make_unique<biconjugate_gradient_stabilised_ocl>(
                    solver_data["tolerance"].get<double>());
            }
            else if (solver_data.find("maximum_iterations") != solver_data.end())
            {
                return std::make_unique<biconjugate_gradient_stabilised_ocl>(
                    solver_data["maximum_iterations"].get<int>());
            }
            return std::make_unique<biconjugate_gradient_stabilised_ocl>();
        }
        else
        {
            if (solver_data.find("tolerance") != solver_data.end()
                && solver_data.find("maximum_iterations") != solver_data.end())
            {
                return std::make_unique<conjugate_gradient_ocl>(solver_data["tolerance"].get<double>(),
                                                                solver_data["maximum_iterations"]
                                                                    .get<int>());
            }
            else if (solver_data.find("tolerance") != solver_data.end())
            {
                return std::make_unique<conjugate_gradient_ocl>(
                    solver_data["tolerance"].get<double>());
            }
            else if (solver_data.find("maximum_iterations") != solver_data.end())
            {
                return std::make_unique<conjugate_gradient_ocl>(
                    solver_data["maximum_iterations"].get<int>());
            }
            return std::make_unique<conjugate_gradient_ocl>();
        }

#else
        throw std::domain_error("iterativeOCL is only available when neon is "
                                "configured with -DENABLE_OCL=1\n");
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
