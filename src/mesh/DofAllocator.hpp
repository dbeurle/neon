
#pragma once

#include "numeric/IndexTypes.hpp"

#include <algorithm>
#include <initializer_list>

namespace neon
{
/**
* Allocates the dof lists from the nodal connectivity vector
* @param nodal_dofs Number of degrees of freedom for each node
* @param nodal_connectivity Vector of nodal coordinates
* @return The global degrees of freedom
*/
[[nodiscard]] inline local_indices allocate_element_dof_list(int32 const nodal_dofs,
                                                             local_indices const& element_nodes)
{
    local_indices element_dof_list;

    element_dof_list.reserve(nodal_dofs * element_nodes.size());

    for (auto const& node : element_nodes)
    {
        for (auto dof_offset = 0; dof_offset < nodal_dofs; ++dof_offset)
        {
            element_dof_list.emplace_back(node * nodal_dofs + dof_offset);
        }
    }
    return element_dof_list;
}

[[nodiscard]] inline std::vector<local_indices> allocate_dof_list(
    int32 const nodal_dofs, std::vector<local_indices> const& element_node_list)
{
    std::vector<local_indices> element_dof_list(element_node_list.size());

    std::transform(std::begin(element_node_list),
                   std::end(element_node_list),
                   std::begin(element_dof_list),
                   [&](auto const& element_nodes) {
                       return allocate_element_dof_list(nodal_dofs, element_nodes);
                   });

    return element_dof_list;
}

/**
   * Allocates the dof lists from the nodal connectivity vector
   * @param nodal_dofs Number of degrees of freedom for each node
   * @param nodal_connectivity Vector of nodal coordinates
   * @return The global degrees of freedom
   */
[[nodiscard]] std::vector<int32> allocate_dof_list(int32 const nodal_dofs,
                                                   std::vector<int32> const& nodal_connectivity);

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
