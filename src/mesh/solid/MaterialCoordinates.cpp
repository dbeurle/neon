/*
 * neon - A finite element solver.
 *
 * For licensing please refer to the LICENSE.md file
 *
 */

#include "MaterialCoordinates.hpp"

namespace neon::solid
{
MaterialCoordinates::MaterialCoordinates(Vector const& initial_coordinates)
    : NodalCoordinates(initial_coordinates), x(initial_coordinates)
{
}

Vector MaterialCoordinates::displacement(List const& local_dofs) const
{
    using namespace ranges;

    Vector localdisp(local_dofs.size());

    for_each(view::zip(view::ints(0), local_dofs), [&](auto const& zip_pair) {
        auto const & [ i, local_dof ] = zip_pair;
        localdisp(i) = x(local_dof) - X(local_dof);
    });
    return localdisp;
}

Matrix MaterialCoordinates::get_configuration(List const& local_nodes, Vector const& configuration) const
{
    auto const lnodes = local_nodes.size();
    Matrix localconf(3, lnodes);

    for (auto lnode = 0; lnode < lnodes; lnode++)
    {
        localconf.col(lnode) = configuration.segment<3>(3 * local_nodes[lnode]);
    }
    return localconf;
}
}
