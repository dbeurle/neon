
#pragma once

#include "mesh/diffusion/Submesh.hpp"

#include "mesh/generic/Dirichlet.hpp"
#include "mesh/diffusion/boundary/SurfaceFlux.hpp"

#include <map>
#include <string>

namespace neon
{
class BasicMesh;

namespace diffusion
{
class femMesh
{
public:
    femMesh(BasicMesh const& basic_mesh,
            Json::Value const& material_data,
            Json::Value const& simulation_data);

    auto active_dofs() const { return material_coordinates->size(); }

    /**
     * Deform the body by updating the displacement x = X + u
     * and update the internal variables with the new deformation and the
     * time step increment
     */
    void update_internal_variables(Vector const& u, double const time_step_size = 0.0);

    /**
     * Update the internal variables if converged, otherwise revert back
     * for next attempted load increment
     */
    void save_internal_variables(bool const have_converged);

    /** Constant access to the sub-meshes */
    std::vector<femSubmesh> const& meshes() const { return submeshes; }

    /** Mutable access to the sub-meshes */
    std::vector<femSubmesh>& meshes() { return submeshes; }

    auto const& dirichlet_boundaries() const { return dirichlet_bcs; }

    auto const& surface_boundaries() const { return surface_bcs; }

    auto const& coordinates() const { return *(material_coordinates.get()); }

protected:
    void check_boundary_conditions(Json::Value const& boundary_data) const;

    void allocate_boundary_conditions(Json::Value const& boundary_data, BasicMesh const& basic_mesh);

    /** Collapse the nodal connectivity arrays from the submesh for a node list */
    List filter_dof_list(std::vector<SubMesh> const& boundary_mesh) const;

protected:
    std::shared_ptr<MaterialCoordinates> material_coordinates;

    std::vector<femSubmesh> submeshes;

    std::map<std::string, std::vector<Dirichlet>> dirichlet_bcs;
    std::map<std::string, std::vector<SurfaceFluxBoundary>> surface_bcs;

    std::unordered_map<std::string, int> const dof_table = {{"x", 0}, {"y", 1}, {"z", 2}};
};
}
}
