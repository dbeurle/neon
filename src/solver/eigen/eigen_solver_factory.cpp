
#include "solver/eigen/eigen_solver_factory.hpp"

#include "io/json.hpp"

#include "solver/eigen/arpack.hpp"
#include "solver/eigen/lanczos_ocl.hpp"
#include "solver/eigen/power_iteration.hpp"

namespace neon
{
std::unique_ptr<eigen_solver> make_eigen_solver(json const& solver_data)
{
    if (solver_data.find("type") == end(solver_data))
    {
        throw std::domain_error("Eigen solver type was not provided.  Please use "
                                "\"power_iteration\", "
                                "\"arpack\" or \"lanczos\"");
    }

    std::int64_t number_of_ev = 10;

    if (solver_data.find("eigenvalues") != end(solver_data))
    {
        number_of_ev = solver_data["eigenvalues"];

        if (number_of_ev <= 0)
        {
            throw std::domain_error("The number of requested eigenvalues must be positive");
        }
    }

    eigen_solver::eigen_spectrum spectrum;

    if (solver_data.find("spectrum") != end(solver_data))
    {
        std::string const& input_spectrum = solver_data["spectrum"];

        if (input_spectrum == "lower")
        {
            spectrum = eigen_solver::eigen_spectrum::lower;
        }
        else if (input_spectrum == "upper")
        {
            spectrum = eigen_solver::eigen_spectrum::upper;
        }
        else
        {
            throw std::domain_error("A spectrum was requested but only \"lower\" and \"upper\" "
                                    "options are supported");
        }
    }

    if (std::string const& type = solver_data["type"]; type == "power_iteration")
    {
        return std::make_unique<power_iteration>(number_of_ev, spectrum);
    }
    else if (type == "lanczos")
    {
        return std::make_unique<lanczos_ocl>(number_of_ev, spectrum);
    }
    else if (type == "arpack")
    {
        return std::make_unique<arpack>(number_of_ev, spectrum);
    }
    return nullptr;
}
}
