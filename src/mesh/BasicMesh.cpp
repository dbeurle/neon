
#include "BasicMesh.hpp"

#include "Exceptions.hpp"

#include "io/json.hpp"

namespace neon
{
BasicMesh::BasicMesh(json const& mesh_file) : NodalCoordinates(mesh_file)
{
    if (!mesh_file.isMember("Elements"))
    {
        throw std::runtime_error("The mesh file is missing the \"Elements\" field");
    }

    for (auto const& mesh : mesh_file["Elements"])
    {
        meshes_map[mesh["Name"].asString()].push_back(mesh);
    }
}

std::vector<Submesh> const& BasicMesh::meshes(std::string const& name) const
{
    auto const found = meshes_map.find(name);

    if (found == meshes_map.end()) throw KeyNotFoundInMap<std::string>(name);

    return found->second;
}
}
