
#include "BasicMesh.hpp"

#include "Exceptions.hpp"

#include "mpi.hpp"

#include <json/value.h>

namespace neon
{
BasicMesh::BasicMesh(Json::Value const& mesh_file) : NodalCoordinates(mesh_file)
{
    if (!mesh_file.isMember("Elements"))
    {
        throw std::runtime_error("The mesh file is missing the \"Elements\" field");
    }

    for (auto const& mesh : mesh_file["Elements"])
    {
        meshes_map[mesh["Name"].asString()].push_back(mesh);
    }

    if (mpi::size() > 1)
    {
        if (!mesh_file.isMember("Interface"))
        {
            throw std::runtime_error("Distributed mesh doesn't have an \"Interface\" field");
        }
        if (!mesh_file.isMember("LocalToGlobalMap"))
        {
            throw std::runtime_error("Distributed mesh doesn't have an \"LocalToGlobalMap\" field");
        }

        // Fill the interface map
        for (auto const& interface : mesh_file["Interface"])
        {
            if (!interface.isMember("Process"))
            {
                throw std::runtime_error("Field \"Process\" required in mesh file");
            }
            if (!interface.isMember("Indices"))
            {
                throw std::runtime_error("Field \"Indices\" required in mesh file");
            }

            auto& interface_nodes = interfaces[interface["Process"].asInt()];

            interface_nodes.reserve(interface["Indices"].size());

            for (auto const& jnode : interface["Indices"])
            {
                interface_nodes.emplace_back(jnode.asInt64());
            }
        }

        nonlocal_indices.reserve(mesh_file["LocalToGlobalMap"].size());

        for (auto const& index : mesh_file["LocalToGlobalMap"])
        {
            nonlocal_indices.emplace_back(index.asInt64());
        }
    }
}

std::vector<Submesh> const& BasicMesh::meshes(std::string const& name) const
{
    if (!this->has(name)) throw KeyNotFoundInMap<std::string>(name);

    return meshes_map.find(name)->second;
}
}
