
#include "Neumann.hpp"

namespace neon
{
Neumann::Neumann(std::vector<List> const& nodal_connectivity,
                 std::vector<List> const& dof_list,
                 std::shared_ptr<MaterialCoordinates>& material_coordinates,
                 Json::Value const& times,
                 Json::Value const& loads)
    : Boundary(times, loads),
      nodal_connectivity(nodal_connectivity),
      dof_list(dof_list),
      material_coordinates(material_coordinates)
{
}
}
