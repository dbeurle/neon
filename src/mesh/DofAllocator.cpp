
#include "DofAllocator.hpp"

#include <range/v3/view.hpp>

namespace neon
{
namespace solid
{
int encode_dof(std::string const& name) { return 0; }
}

std::vector<List> allocate_dof_list(int const nodal_dofs, std::vector<List> const& nodal_connectivity)
{
    using namespace ranges;
    return nodal_connectivity | view::transform([=](auto const& node_list) {
               return view::for_each(node_list, [=](int const local_node) {
                   return view::ints(0, nodal_dofs) | view::transform([=](int const nodal_dof) {
                              return local_node * nodal_dofs + nodal_dof;
                          });
               });
           });
}

std::vector<List> filter_dof_list(int const nodal_dofs,
                                  int const dof_offset,
                                  std::vector<List> const& nodal_connectivity)
{
    using namespace ranges;
    return nodal_connectivity | view::transform([=](auto const& node_list) {
               return node_list
                      | view::transform([=](auto const& n) { return n * nodal_dofs + dof_offset; });
           });
}
}
