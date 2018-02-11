
#pragma once

#include "numeric/index_types.hpp"

namespace neon
{
/**
 * Allocates the dof lists from the nodal connectivity vector
 * @param nodal_dofs Number of degrees of freedom for each node
 * @param nodal_connectivity vector of nodal coordinates
 * @return The global degrees of freedom
 */
[[nodiscard]] std::vector<local_indices> allocate_dof_list(
    std::int32_t const nodal_dofs, std::vector<local_indices> const& element_node_list);

/**
 * Allocates the dof lists from the nodal connectivity vector
 * @param nodal_dofs Number of degrees of freedom for each node
 * @param nodal_connectivity vector of nodal coordinates
 * @return The global degrees of freedom
 */

[[nodiscard]] local_indices allocate_element_dofs(std::int32_t const nodal_dofs,
                                                  local_indices const& element_nodes);

/**
 * This function accepts the nodal connectivity of the mesh, multiplies each
 * entry by the nodal_dofs and adds the offset.  This function is intended for
 * use with the boundary classes, where each boundary class holds a dof_list
 * with the dofs associated only with the particular dof.
 */
[[nodiscard]] std::vector<local_indices> filter_dof_list(
    int const nodal_dofs, int const dof_offset, std::vector<local_indices> const& nodal_connectivity);
}
