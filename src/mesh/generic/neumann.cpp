
#include "neumann.hpp"

namespace neon
{
neumann::neumann(indices nodal_connectivity,
                 indices dof_list,
                 std::shared_ptr<material_coordinates>& coordinates,
                 json const& times,
                 json const& loads)
    : vector_contribution{times, loads},
      nodal_connectivity{nodal_connectivity},
      dof_list{dof_list},
      coordinates{coordinates}
{
}

neumann::neumann(indices nodal_connectivity,
                 indices dof_list,
                 std::shared_ptr<material_coordinates>& coordinates,
                 json const& boundary,
                 std::string const& name,
                 double const generate_time_step)
    : vector_contribution{boundary, name, generate_time_step},
      nodal_connectivity{nodal_connectivity},
      dof_list{dof_list},
      coordinates{coordinates}
{
}
}
