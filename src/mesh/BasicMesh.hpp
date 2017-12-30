
#pragma once

#include "mesh/NodalCoordinates.hpp"
#include "mesh/Submesh.hpp"

#include <map>

namespace neon
{
/**
 * BasicMesh is a basic mesh definition, with nodes and associated volume
 * meshes. This container does not define a notion of boundary or volume meshes.
 * It holds nodes and element type collections, defining type and nodal
 * connectivities.
 *
 * This class models the following mesh structure:
   \verbatim
                         parent mesh
                              |
              name0 ----------------------- name1
                |                             |
            mesh_type2 ---- mesh_type1
                |               |
            hexahedrons      prisms
   \endverbatim
 *
 */
class BasicMesh : public NodalCoordinates
{
public:
    BasicMesh(Json::Value const& mesh_file);

    /** @return mesh matching a specific name */
    [[nodiscard]] std::vector<Submesh> const& meshes(std::string const& name) const;

    /** @return true if the name exists in this store */
    bool has(std::string const& name) const { return meshes_map.find(name) != meshes_map.end(); }

    /** @return the process and the shared global node indices */
    auto process_interfaces() const { return interfaces; }

    /** @return local mapping  */
    auto local_to_nonlocal_ordering() const { return nonlocal_indices; }

protected:
    std::map<std::string, std::vector<Submesh>> meshes_map;

    std::unordered_map<std::int32_t, std::vector<int64>> interfaces;

    std::vector<int64> nonlocal_indices;
};
}
