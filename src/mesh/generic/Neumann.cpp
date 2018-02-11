
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
}
