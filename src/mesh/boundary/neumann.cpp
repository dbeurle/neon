
#include "neumann.hpp"

#include <utility>

namespace neon
{
neumann::neumann(indices node_indices,
                 indices dof_indices,
                 std::shared_ptr<material_coordinates>& coordinates,
                 json const& times,
                 json const& loads)
    : vector_contribution{times, loads},
      node_indices{std::move(node_indices)},
      dof_indices{std::move(dof_indices)},
      coordinates{coordinates}
{
}

neumann::neumann(indices node_indices,
                 indices dof_indices,
                 std::shared_ptr<material_coordinates>& coordinates,
                 json const& boundary,
                 std::string const& name,
                 double const generate_time_step)
    : vector_contribution{boundary, name, generate_time_step},
      node_indices{std::move(node_indices)},
      dof_indices{std::move(dof_indices)},
      coordinates{coordinates}
{
}
}
