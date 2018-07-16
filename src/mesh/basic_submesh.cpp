
#include "basic_submesh.hpp"

#include "mesh/node_ordering_adapter.hpp"
#include "exceptions.hpp"
#include "io/json.hpp"

#include <set>

namespace neon
{
basic_submesh::basic_submesh(json const& mesh)
{
    // Error checking for empty fields
    if (!mesh.count("Name"))
    {
        throw std::domain_error("The element group in the mesh file is missing the "
                                "\"Name\" field");
    }
    if (!mesh.count("Type"))
    {
        throw std::domain_error("The element group in the mesh file is missing the "
                                "\"Type\" field");
    }
    if (!mesh.count("NodalConnectivity"))
    {
        throw std::domain_error("The element group in the mesh file is missing the "
                                "\"NodalConnectivity\" field");
    }
    if (mesh["NodalConnectivity"].size() == 0 || mesh["NodalConnectivity"][0].size() == 0)
    {
        throw std::domain_error("The element group in the mesh file is empty");
    }

    m_topology = gmsh_type_to_enum(mesh["Type"]);

    node_indices.resize(mesh["NodalConnectivity"][0].size(), mesh["NodalConnectivity"].size());

    {
        std::int64_t element_counter{0};

        for (auto const& mesh_connectivity : mesh["NodalConnectivity"])
        {
            std::int64_t node_counter{0};
            for (auto const& node : mesh_connectivity)
            {
                node_indices(node_counter++, element_counter) = node;
            }
            ++element_counter;
        }
    }
    convert_from_gmsh(node_indices, m_topology);
}

std::vector<std::int32_t> basic_submesh::unique_node_indices() const
{
    std::set<std::int32_t> unique_set;

    std::copy_n(node_indices.data(), node_indices.size(), std::inserter(unique_set, end(unique_set)));

    return {begin(unique_set), end(unique_set)};
}
}
