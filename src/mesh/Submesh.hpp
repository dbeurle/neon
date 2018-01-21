
#pragma once

#include "mesh/ElementTopology.hpp"
#include "numeric/IndexTypes.hpp"

#include <json/forwards.h>

namespace neon
{
/** Submesh stores connectivities for one particular element group */
class Submesh
{
public:
    /** Construct using a json object */
    Submesh(json const& mesh);

    /** @return the number of elements */
    auto elements() const { return nodal_connectivity.size(); }

    /** @return The element topology for this mesh */
    ElementTopology const& topology() const { return element_topology; }

    /** @return a list of the element nodal connectivities */
    List const& local_node_list(int const element) const { return nodal_connectivity[element]; }

    /** @return the number of nodes in the element */
    auto nodes_per_element() const { return nodal_connectivity.at(0).size(); }

    /** @return a list of unique nodal connectivities */
    List unique_connectivities() const;

    /** @return a vector of lists, with each list giving the element nodes */
    auto const& connectivities() const { return nodal_connectivity; }

protected:
    ElementTopology element_topology;
    std::vector<List> nodal_connectivity;
};
}
