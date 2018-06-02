
#pragma once

#include "linear_solver.hpp"
#include "io/json_forward.hpp"

#include <memory>

namespace neon
{
std::unique_ptr<linear_solver> make_linear_solver(json const& solver_data,
                                                  bool const is_symmetric = true);
}
