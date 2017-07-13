
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

    /** Reset the boundary conditions */
    void internal_restart(Json::Value const& simulation_data);

    /**
     * Deform the body by updating the displacement x = X + u
     * and update the internal variables with the new deformation and the
     * time step increment
     */
    void update_internal_variables(Vector const& u, double const Δt = 1.0);

    /**
     * Update the internal variables if converged, otherwise revert back
     * for next attempted load increment
     */
    void save_internal_variables(bool const have_converged);

    /** Constant access to the sub-meshes */
    std::vector<femSubmesh> const& meshes() const { return submeshes; }

    /** Mutable access to the sub-meshes */
    std::vector<femSubmesh>& meshes() { return submeshes; }

    auto const& dirichlet_boundary_map() const { return dirichlet_boundaries; }

    auto const& coordinates() const { return *(material_coordinates.get()); }

    void write(int filename_append = -1) const;

protected:
    void check_boundary_conditions(Json::Value const& boundary_data) const;

    void allocate_boundary_conditions(Json::Value const& boundary_data, BasicMesh const& basic_mesh);

    /** \sa internal_restart */
    void reallocate_boundary_conditions(Json::Value const& boundary_data);

    /** Collapse the nodal connectivity arrays from the submesh for a node list */
    List filter_dof_list(std::vector<SubMesh> const& boundary_mesh) const;

protected:
    std::shared_ptr<MaterialCoordinates> material_coordinates;
    std::vector<femSubmesh> submeshes;

    std::map<std::string, std::vector<Dirichlet>> dirichlet_boundaries;

    const std::unordered_map<std::string, int> dof_table = {{"x", 0}, {"y", 1}, {"z", 2}};
};
}
}
