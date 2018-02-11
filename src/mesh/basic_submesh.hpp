
#pragma once

#include "mesh/element_topology.hpp"
#include "numeric/index_types.hpp"

#include "io/json_forward.hpp"

namespace neon
{
/** basic_submesh stores connectivities for one particular element group */
class basic_submesh
{
public:
    /** Construct using a json object */
    basic_submesh(json const& mesh);

    /** @return the number of elements */
    auto elements() const { return nodal_connectivity.size(); }

    /** @return The element topology for this mesh */
    element_topology const& topology() const { return m_topology; }

    /** @return a list of the element nodal connectivities */
    local_indices const& local_node_list(int const element) const { return nodal_connectivity[element]; }

    /** @return the number of nodes in the element */
    auto nodes_per_element() const { return nodal_connectivity.at(0).size(); }

    /** @return a list of unique nodal connectivities */
    local_indices unique_connectivities() const;

    /** @return a vector of lists, with each list giving the element nodes */
    auto const& connectivities() const { return nodal_connectivity; }

protected:
    element_topology m_topology;
    std::vector<local_indices> nodal_connectivity;
};
}
