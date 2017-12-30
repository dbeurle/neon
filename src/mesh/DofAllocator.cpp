
#include "DofAllocator.hpp"

#include <range/v3/view/for_each.hpp>
#include <range/v3/view/iota.hpp>
#include <range/v3/view/transform.hpp>

#include <algorithm>

namespace neon
{
std::vector<int64> allocate_element_dofs(int32 const nodal_dofs,
                                         std::vector<int64> const& element_nodes)
{
    std::vector<int64> element_dof_list;

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

std::vector<std::vector<int64>> allocate_dof_list(
    int32 const nodal_dofs, std::vector<std::vector<int64>> const& element_node_list)
{
    std::vector<std::vector<int64>> element_dof_list(element_node_list.size());

    std::transform(std::begin(element_node_list),
                   std::end(element_node_list),
                   std::begin(element_dof_list),
                   [&](auto const& element_nodes) {
                       return allocate_element_dofs(nodal_dofs, element_nodes);
                   });

    return element_dof_list;
}

std::vector<std::vector<int64>> filter_dof_list(int32 const nodal_dofs,
                                                int32 const dof_offset,
                                                std::vector<std::vector<int64>> const& nodal_connectivity)
{
    using namespace ranges;
    return nodal_connectivity | view::transform([=](auto const& node_list) {
               return node_list
                      | view::transform([=](auto const& n) { return n * nodal_dofs + dof_offset; });
           });
}
}
