
#pragma once

#include "LinearSolver.hpp"

#include <memory>

#include <json/forwards.h>

namespace neon
{
std::unique_ptr<LinearSolver> make_linear_solver(json const& solver_data,
                                                 bool const is_symmetric = true);
}
