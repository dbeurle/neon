
#pragma once

#include "numeric/DenseTypes.hpp"

#include <range/v3/view.hpp>

namespace neon
{
/**
 * Allocates the dof lists from the nodal connectivity vector
 * @param nodal_dofs Number of degrees of freedom for each node
 * @param nodal_connectivity Vector of nodal coordinates
 * @return The global degrees of freedom
 */
std::vector<List> allocate_dof_list(int const nodal_dofs,
                                    std::vector<List> const& nodal_connectivity)
{
    using namespace ranges;
    return nodal_connectivity | view::transform([=](auto const& node_list) {
               return view::for_each(node_list, [=](int const local_node) {
                   return view::ints(0, nodal_dofs)
                          | view::transform([=](int const nodal_dof) {
                                return local_node * nodal_dofs + nodal_dof;
                            });
               });
           });
}

/**
 * This function accepts the nodal connectivity of the mesh, multiplies each
 * entry by the nodal_dofs and adds the offset.  This function is intended for
 * use with the boundary classes, where each boundary class holds a dof_list
 * with the dofs associated only with it's
 */
std::vector<List> filter_dof_list(int const nodal_dofs,
                                  int const dof_offset,
                                  std::vector<List> const& nodal_connectivity)
{
    using namespace ranges;
    return nodal_connectivity | view::transform([=](auto const& node_list) {
               return node_list | view::transform([=](auto const& n) {
                          return n * nodal_dofs + dof_offset;
                      });
           });
}
}
