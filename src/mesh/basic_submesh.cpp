
#include "basic_submesh.hpp"

#include "mesh/node_ordering_adapter.hpp"

#include "Exceptions.hpp"

#include "io/json.hpp"

#include <range/v3/action/join.hpp>
#include <range/v3/action/sort.hpp>
#include <range/v3/action/unique.hpp>

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

    connectivity.resize(mesh["NodalConnectivity"][0].size(), mesh["NodalConnectivity"].size());

    {
        std::int64_t element_counter{0};

        for (auto const& mesh_connectivity : mesh["NodalConnectivity"])
        {
            std::int64_t node_counter{0};
            for (auto const& node : mesh_connectivity)
            {
                connectivity(node_counter++, element_counter) = node;
            }
            ++element_counter;
        }
    }
    convert_from_gmsh(connectivity, m_topology);
}

std::vector<std::int32_t> basic_submesh::unique_connectivity() const
{
    std::set<std::int32_t> unique_set;

    std::copy_n(connectivity.data(),
                connectivity.size(),
                std::inserter(unique_set, std::end(unique_set)));

    return {std::begin(unique_set), std::end(unique_set)};
}
}
