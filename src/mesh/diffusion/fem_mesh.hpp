
#pragma once

#include "mesh/diffusion/fem_submesh.hpp"

#include "mesh/diffusion/boundary/surface_boundary.hpp"
#include "mesh/generic/dirichlet.hpp"

#include <map>
#include <string>

namespace neon
{
class basic_mesh;

namespace diffusion
{
class fem_mesh
{
public:
    fem_mesh(basic_mesh const& basic_mesh, json const& material_data, json const& mesh_data);

    auto active_dofs() const { return mesh_coordinates->size(); }

    /**
     * Deform the body by updating the displacement x = X + u
     * and update the internal variables with the new deformation and the
     * time step increment
     */
    void update_internal_variables(vector const& u, double const time_step_size = 0.0);

    /**
     * Update the internal variables if converged, otherwise revert back
     * for next attempted load increment
     */
    void save_internal_variables(bool const have_converged);

    /** Constant access to the sub-meshes */
    std::vector<fem_submesh> const& meshes() const { return submeshes; }

    /** Mutable access to the sub-meshes */
    std::vector<fem_submesh>& meshes() { return submeshes; }

    auto const& dirichlet_boundaries() const { return dirichlet_bcs; }

    auto const& surface_boundaries() const { return boundary_meshes; }

    [[deprecated]] auto const& coordinates() const { return *mesh_coordinates; }

    auto const& geometry() const { return *mesh_coordinates; }

protected:
    void check_boundary_conditions(json const& boundary_data) const;

    void allocate_boundary_conditions(json const& boundary_data, basic_mesh const& basic_mesh);

protected:
    std::shared_ptr<material_coordinates> mesh_coordinates;

    std::vector<fem_submesh> submeshes;

    std::map<std::string, std::vector<dirichlet>> dirichlet_bcs;
    std::map<std::string, std::vector<boundary_mesh>> boundary_meshes;

    std::unordered_map<std::string, int> const dof_table = {{"x", 0}, {"y", 1}, {"z", 2}};
};
}
}
