
#pragma once

#include "mesh/NodalCoordinates.hpp"
#include "mesh/SubMesh.hpp"

#include <map>

namespace neon
{
/**
 * BasicMesh is a basic mesh definition, with nodes and associated volume
 * meshes. This container does not define a notion of boundary or volume meshes.
 * It holds nodes and element type collections, defining type and nodal
 * connectivities.
 */
class BasicMesh : public NodalCoordinates
{
public:
    BasicMesh(Json::Value const& mesh_file);

    /** @return mesh matching a specific name */
    std::vector<SubMesh> const& meshes(std::string const& name) const;

protected:
    std::map<std::string, std::vector<SubMesh>> meshes_map;
};
}
