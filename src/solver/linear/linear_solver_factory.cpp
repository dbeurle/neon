
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
#include <set>

namespace neon
{
std::unique_ptr<linear_solver> make_linear_solver(json const& solver_data, bool const is_symmetric)
{
    if (solver_data.find("type") == end(solver_data))
    {
        throw std::domain_error("A \"linear_solver\" \"type\" must be specified");
    }

    std::string const& solver_name = solver_data["type"];

    {
        std::set<std::string> names{"PaStiX", "MUMPS", "direct", "iterative"};

        if (!std::binary_search(begin(names), end(names), solver_name))
        {
            throw std::domain_error("Linear solver " + solver_name
                                    + " is not recognised.  Please use \"PaStiX\", "
                                      "\"MUMPS\", \"direct\" or \"iterative\"");
        }
    }

    // If CUDA is enabled we need to check if there are valid compute
    // devices available at runtime.  If OpenCL is available, then we prefer this
    // and check for valid compute devices.  Otherwise we will fall back to a
    // multithreaded CPU implementation.
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
    else if (solver_name == "direct")
    {
        if (is_symmetric) return std::make_unique<SparseLLT>();

        return std::make_unique<SparseLU>();
    }
    else if (solver_name == "iterative")
    {
        // Alias the possible solvers
#if defined ENABLE_CUDA
        using cg_type = conjugate_gradient_cuda;
#elif defined ENABLE_OCL
        using cg_type = conjugate_gradient_ocl;
#else
        using cg_type = conjugate_gradient;
#endif

#if defined ENABLE_CUDA
        using bicg_type = biconjugate_gradient_stabilised_cuda;
#elif defined ENABLE_OCL
        using bicg_type = biconjugate_gradient_stabilised_ocl;
#else
        using bicg_type = biconjugate_gradient_stabilised;
#endif

        if (solver_data.find("tolerance") != solver_data.end()
            && solver_data.find("maximum_iterations") != solver_data.end())
        {
            double const tolerance = solver_data["tolerance"];
            std::int32_t const maximum_iterations = solver_data["maximum_iterations"];

            return is_symmetric ? std::make_unique<cg_type>(tolerance, maximum_iterations)
                                : std::make_unique<bicg_type>(tolerance, maximum_iterations);
        }
        else if (solver_data.find("tolerance") != solver_data.end())
        {
            double const tolerance = solver_data["tolerance"];

            return is_symmetric ? std::make_unique<cg_type>(tolerance)
                                : std::make_unique<bicg_type>(tolerance);
        }
        else if (solver_data.find("maximum_iterations") != solver_data.end())
        {
            std::int32_t const maximum_iterations = solver_data["maximum_iterations"];

            return is_symmetric ? std::make_unique<cg_type>(maximum_iterations)
                                : std::make_unique<bicg_type>(maximum_iterations);
        }
        return is_symmetric ? std::make_unique<cg_type>() : std::make_unique<bicg_type>();
    }
    else
    {
        throw std::domain_error("Did not find a linear solver type.  Did you try specifying "
                                "\"type\" for the linear solver name?\n");
    }
    return nullptr;
}
}
