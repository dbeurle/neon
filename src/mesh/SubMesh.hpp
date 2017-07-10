
#pragma once

#include <json/forwards.h>

#include "mesh/ElementTopology.hpp"
#include "mesh/NodeOrderingAdapter.hpp"
#include "numeric/DenseTypes.hpp"

namespace neon
{
class NodalCoordinates;

/** SubMesh stores connectivities for one particular element group */
class SubMesh
{
public:
    /** Construct using a json object */
    SubMesh(Json::Value const& mesh);

    /** @return the number of elements */
    auto elements() const { return nodal_connectivity.size(); }

    /** @return The element topology for this mesh */
    ElementTopology topology() const { return element_topology; }

    /** @return a list of the element nodal connectivities */
    List const& local_node_list(int const element) const { return nodal_connectivity[element]; }

    /** @return the number of nodes in the element */
    auto nodes_per_element() const { return nodal_connectivity.at(0).size(); }

    /** @return a list of unique nodal connectivities */
    List unique_connectivities() const;

    /** @return a vector of lists, with each list giving the element nodes */
    auto const& connectivities() const { return nodal_connectivity; }

protected:
    /** @return coordinates given a Nodes class and  */
    Matrix gather_coordinates(NodalCoordinates const& nodal_coordinates, int const element) const;

protected:
    ElementTopology element_topology;
    NodeOrderingAdapter adapter;
    std::vector<List> nodal_connectivity;
};
}
