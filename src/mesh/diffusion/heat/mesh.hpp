
#pragma once

/// @file

#include "mesh/diffusion/heat/submesh.hpp"

#include "mesh/diffusion/heat/boundary/surface_boundary.hpp"
#include "mesh/boundary/dirichlet.hpp"
#include "io/file_output.hpp"

#include <map>
#include <memory>

namespace neon
{
class basic_mesh;

namespace diffusion
{
class mesh
{
public:
    mesh(basic_mesh const& basic_mesh, json const& material_data, json const& mesh_data);

    auto active_dofs() const { return coordinates->size(); }

    /// Deform the body by updating the displacement x = X + u
    /// and update the internal variables with the new deformation and the
    /// time step increment
    void update_internal_variables(vector const& u, double const time_step_size = 0.0);

    /// Update the internal variables if converged, otherwise revert back
    /// for next attempted load increment
    void save_internal_variables(bool const have_converged);

    /// Constant access to the sub-meshes
    std::vector<submesh> const& meshes() const { return submeshes; }

    /// Mutable access to the sub-meshes
    std::vector<submesh>& meshes() { return submeshes; }

    auto const& dirichlet_boundaries() const { return dirichlet_bcs; }

    auto const& surface_boundaries() const { return boundary_meshes; }

    auto const& geometry() const { return *coordinates; }

    /// Write out results to file
    void write(std::int32_t const time_step, double const current_time);

protected:
    void check_boundary_conditions(json const& boundary_data) const;

    void allocate_variable_names();

    void allocate_boundary_conditions(json const& boundary_data, basic_mesh const& basic_mesh);

protected:
    std::shared_ptr<material_coordinates> coordinates;

    std::vector<submesh> submeshes;

    std::map<std::string, std::vector<dirichlet>> dirichlet_bcs;
    std::map<std::string, std::vector<boundary_mesh>> boundary_meshes;

    std::unordered_map<std::string, int> const dof_table = {{"x", 0}, {"y", 1}, {"z", 2}};

    /// File output handle
    std::unique_ptr<io::file_output> writer;

    /// Output variables
    std::vector<variable::types> output_variables;

    vector temperature;
};
}
}
