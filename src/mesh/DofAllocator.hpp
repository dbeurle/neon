
#pragma once

#include "numeric/IndexTypes.hpp"


namespace neon
{
/**
 * Allocates the dof lists from the nodal connectivity vector
 * @param nodal_dofs Number of degrees of freedom for each node
 * @param nodal_connectivity Vector of nodal coordinates
 * @return The global degrees of freedom
 */
[[nodiscard]] std::vector<List> allocate_dof_list(int const nodal_dofs,
                                                  std::vector<List> const& nodal_connectivity);

/**
 * This function accepts the nodal connectivity of the mesh, multiplies each
 * entry by the nodal_dofs and adds the offset.  This function is intended for
 * use with the boundary classes, where each boundary class holds a dof_list
 * with the dofs associated only with the particular dof.
 */
[[nodiscard]] std::vector<List> filter_dof_list(int const nodal_dofs,
                                                int const dof_offset,
                                                std::vector<List> const& nodal_connectivity);
}
