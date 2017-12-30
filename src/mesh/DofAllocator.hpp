
#pragma once

#include "numeric/IndexTypes.hpp"

namespace neon
{
/**
 * Allocates the dof lists from the nodal connectivity vector
 * @param nodal_dofs Number of degrees of freedom for each node
 * @param nodal_connectivity vector of nodal coordinates
 * @return The global degrees of freedom
 */
[[nodiscard]] std::vector<std::vector<int64>> allocate_dof_list(
    int32 const nodal_dofs, std::vector<std::vector<int64>> const& element_node_list);

/**
 * Allocates the dof lists from the nodal connectivity vector
 * @param nodal_dofs Number of degrees of freedom for each node
 * @param nodal_connectivity vector of nodal coordinates
 * @return The global degrees of freedom
 */

[[nodiscard]] std::vector<int64> allocate_element_dofs(int32 const nodal_dofs,
                                                       std::vector<int64> const& element_nodes);

/**
 * This function accepts the nodal connectivity of the mesh, multiplies each
 * entry by the nodal_dofs and adds the offset.  This function is intended for
 * use with the boundary classes, where each boundary class holds a dof_list
 * with the dofs associated only with the particular dof.
 */
[[nodiscard]] std::vector<std::vector<int64>> filter_dof_list(
    int32 const nodal_dofs,
    int32 const dof_offset,
    std::vector<std::vector<int64>> const& nodal_connectivity);
}
