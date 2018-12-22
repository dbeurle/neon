
#pragma once

/// @file

#include "mesh/nodal_coordinates.hpp"
#include "mesh/basic_submesh.hpp"

#include <map>

namespace neon
{
/// basic_mesh is a basic mesh definition, with nodes and associated volume
/// meshes. This container does not define a notion of boundary or volume meshes.
/// It holds nodes and element type collections, defining type and nodal
/// connectivities.
class basic_mesh : public nodal_coordinates
{
public:
    basic_mesh(json const& mesh_file);

    /// \return mesh matching a specific name
    [[nodiscard]] std::vector<basic_submesh> const& meshes(std::string const& name) const;

protected:
    std::map<std::string, std::vector<basic_submesh>> meshes_map;
};
}
