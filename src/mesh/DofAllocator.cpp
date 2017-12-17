
#include "DofAllocator.hpp"

#include <range/v3/view/for_each.hpp>
#include <range/v3/view/iota.hpp>
#include <range/v3/view/transform.hpp>

namespace neon
{
std::vector<local_indices> filter_dof_list(int32 const nodal_dofs,
                                           int32 const dof_offset,
                                           std::vector<local_indices> const& nodal_connectivity)
{
    using namespace ranges;
    return nodal_connectivity | view::transform([=](auto const& node_list) {
               return node_list
                      | view::transform([=](auto const& n) { return n * nodal_dofs + dof_offset; });
           });
}
}
