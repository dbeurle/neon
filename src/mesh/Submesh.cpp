
#include "Submesh.hpp"

#include "Exceptions.hpp"

#include <json/value.h>

#include <range/v3/action/join.hpp>
#include <range/v3/action/sort.hpp>
#include <range/v3/action/unique.hpp>

namespace neon
{
Submesh::Submesh(Json::Value const& mesh)
{
    // Error checking for empty fields
    if (!mesh.isMember("Name"))
    {
        throw std::runtime_error("The element group in the mesh file is missing the "
                                 "\"Name\" field");
    }
    if (!mesh.isMember("Type"))
    {
        throw std::runtime_error("The element group in the mesh file is missing the "
                                 "\"Type\" field");
    }
    if (!mesh.isMember("NodalConnectivity"))
    {
        throw std::runtime_error("The element group in the mesh file is missing the "
                                 "\"NodalConnectivity\" field");
    }

    element_topology = adapter.gmsh_type_to_enum(mesh["Type"].asInt());

    nodal_connectivity.reserve(mesh["NodalConnectivity"].size());

    for (auto const& mesh_connectivity : mesh["NodalConnectivity"])
    {
        nodal_connectivity.push_back(List());
        nodal_connectivity.back().reserve(mesh_connectivity.size());

        for (auto const& node : mesh_connectivity)
        {
            nodal_connectivity.back().push_back(node.asInt());
        }
    }
    adapter.convert_from_gmsh(nodal_connectivity, element_topology);
}

List Submesh::unique_connectivities() const
{
    using namespace ranges;
    return std::ref(nodal_connectivity) | action::join | action::sort | action::unique;
}
}
