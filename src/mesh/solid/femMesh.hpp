
#pragma once

#include <json/forwards.h>

#include "numeric/DenseTypes.hpp"

#include "mesh/solid/MaterialCoordinates.hpp"
#include "mesh/solid/Submesh.hpp"

#include "mesh/solid/boundary/Dirichlet.hpp"

namespace neon
{
class BasicMesh;

namespace solid
{
class femMesh
{
public:
    femMesh(BasicMesh const& basic_mesh,
            Json::Value const& material_data,
            Json::Value const& simulation_data);

    int active_dofs() const;

    /**
     * Deform the body by updating the displacement x = X + u
     * and update the internal variables with the new deformation
     */
    void update_internal_variables(Vector const& u);

    /** Constant access to the sub-meshes */
    std::vector<femSubmesh> const& meshes() const { return submeshes; }

    /** Mutable access to the sub-meshes */
    std::vector<femSubmesh>& meshes() { return submeshes; }

    auto const& dirichlet_boundary_map() const { return dirichlet_boundaries; }

    auto const& coordinates() const { return *(material_coordinates.get()); }

    void write(int filename_append = -1) const;

protected:
    std::shared_ptr<MaterialCoordinates> material_coordinates;
    std::vector<femSubmesh> submeshes;

    std::map<std::string, std::vector<Dirichlet>> dirichlet_boundaries;

    const std::unordered_map<std::string, int> dof_table = {{"x", 0}, {"y", 1}, {"z", 2}};
};
}
}
