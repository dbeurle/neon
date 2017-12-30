
#include "neumann.hpp"

namespace neon::boundary
{
neumann::neumann(std::vector<std::vector<int64>> const& nodal_connectivity,
                 std::vector<std::vector<int64>> const& dof_list,
                 std::shared_ptr<MaterialCoordinates>& material_coordinates,
                 Json::Value const& times,
                 Json::Value const& loads)
    : vector_contribution(times, loads),
      nodal_connectivity(nodal_connectivity),
      dof_list(dof_list),
      material_coordinates(material_coordinates)
{
}
}
