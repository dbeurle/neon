
#include "Neumann.hpp"

namespace neon
{
Neumann::Neumann(std::vector<local_indices> const& nodal_connectivity,
                 std::vector<local_indices> const& dof_list,
                 std::shared_ptr<material_coordinates>& coordinates,
                 json const& times,
                 json const& loads)
    : VectorContribution(times, loads),
      nodal_connectivity(nodal_connectivity),
      dof_list(dof_list),
      coordinates(coordinates)
{
}

Neumann::Neumann(std::vector<local_indices> const& nodal_connectivity,
                 std::vector<local_indices> const& dof_list,
                 std::shared_ptr<material_coordinates>& coordinates,
                 json const& boundary,
                 std::string const& name,
                 double const generate_time_step)
    : VectorContribution(boundary, name, generate_time_step),
      nodal_connectivity(nodal_connectivity),
      dof_list(dof_list),
      coordinates(coordinates)
{
}
}
