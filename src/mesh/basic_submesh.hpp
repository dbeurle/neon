
#pragma once

#include "mesh/element_topology.hpp"
#include "numeric/index_types.hpp"

#include "io/json_forward.hpp"

namespace neon
{
/// basic_submesh stores nodal indices and element typology for an element group
class basic_submesh
{
public:
    /// Construct using a json object
    basic_submesh(json const& mesh);

    /// \return the number of elements
    auto elements() const { return node_indices.cols(); }

    /// \return the number of nodes in the element
    auto nodes_per_element() const { return node_indices.rows(); }

    /// \return The element topology for this mesh
    element_topology const& topology() const { return m_topology; }

    /// \return a list of the element nodal connectivities
    index_view local_node_view(std::int64_t const element) const
    {
        return node_indices(Eigen::placeholders::all, element);
    }

    /// \return a list of unique nodal connectivities
    std::vector<int32_t> unique_node_indices() const;

    /// \return a two dimensional array with element nodes
    auto const& all_node_indices() const { return node_indices; }

protected:
    element_topology m_topology;

    indices node_indices; /// A column is one element, rows are the element
};
}
