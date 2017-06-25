
#include "BasicMesh.hpp"

#include "PreprocessorExceptions.hpp"

#include <json/json.h>

namespace neon
{
BasicMesh::BasicMesh(Json::Value const& mesh_file) : NodalCoordinates(mesh_file)
{
    if (mesh_file["Elements"].empty()) throw EmptyFieldException("Elements");

    for (auto const& mesh : mesh_file["Elements"])
    {
        meshes_map[mesh["Name"].asString()].push_back(mesh);
    }
}

std::vector<SubMesh> const& BasicMesh::meshes(std::string const& name) const
{
    auto const found = meshes_map.find(name);

    if (found == meshes_map.end()) throw KeyNotFoundInMap<std::string>(name);

    return found->second;
}
}
