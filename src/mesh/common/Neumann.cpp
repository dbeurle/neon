
#include "Neumann.hpp"

namespace neon
{
Neumann::Neumann(std::vector<List> const& nodal_connectivity,
                 std::vector<List> const& dof_list,
                 double const prescribed_value,
                 bool const is_load_ramped)
    : Boundary(prescribed_value, is_load_ramped),
      nodal_connectivity(nodal_connectivity),
      dof_list(dof_list)
{
}
}
