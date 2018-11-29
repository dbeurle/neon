
#pragma once

#include "io/json_forward.hpp"
#include "solver/eigen/eigen_solver.hpp"

namespace neon
{
std::unique_ptr<eigen_solver> make_eigen_solver(json const& solver_data);
}
