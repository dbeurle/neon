
#pragma once

#include "LinearSolver.hpp"

#include <memory>

#include "io/json_forward.hpp"

namespace neon
{
std::unique_ptr<LinearSolver> make_linear_solver(json const& solver_data,
                                                 bool const is_symmetric = true);
}
