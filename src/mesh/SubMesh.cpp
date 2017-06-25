
#include "SubMesh.hpp"

#include "NodalCoordinates.hpp"
#include "NodeOrderingAdapter.hpp"
#include "PreprocessorExceptions.hpp"

#include <json/json.h>
#include <range/v3/action.hpp>

namespace neon
{
SubMesh::SubMesh(Json::Value const& mesh)
{
    // Error checking for empty fields
    if (mesh["Name"].empty()) throw EmptyFieldException("Name");
    if (mesh["Type"].empty()) throw EmptyFieldException("Type");
    if (mesh["NodalConnectivity"].empty()) throw EmptyFieldException("NodalConnectivity");

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

List SubMesh::unique_connectivities() const
{
    using namespace ranges;
    return std::ref(nodal_connectivity) | action::join | action::sort | action::unique;
}

Matrix SubMesh::gather_coordinates(NodalCoordinates const& nodal_coordinates, int const element) const
{
    return nodal_coordinates[local_node_list(element)];
}
}
