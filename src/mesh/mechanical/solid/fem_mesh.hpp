
#pragma once

#include "mesh/generic/dirichlet.hpp"
#include "mesh/mechanical/solid/boundary/nonfollower_load.hpp"
#include "mesh/mechanical/solid/fem_submesh.hpp"
#include "io/file_output.hpp"

#include <map>

namespace neon
{
class basic_mesh;
}

namespace neon::mechanical::solid
{
class fem_mesh
{
public:
    using internal_variable_type = fem_submesh::internal_variable_type;

    /// Alias traits to submesh
    using traits = fem_submesh::traits;

public:
    fem_mesh(basic_mesh const& basic_mesh,
             json const& material_data,
             json const& simulation_data,
             double const generate_time_step);

    /// The number of active degrees of freedom in this mesh
    [[nodiscard]] auto active_dofs() const { return traits::dofs_per_node * coordinates->size(); }

    /// Checks the boundary conditions and constitutive model to ensure
    /// resulting matrix from this mesh is symmetric.  \sa LinearSolver
    [[nodiscard]] bool is_symmetric() const;

    /// Update the internal forces for printing out reaction forces
    void update_internal_forces(vector const& fint) { reaction_forces = -fint; }

    /// \return nodal reaction forces
    vector const& nodal_reaction_forces() const { return reaction_forces; }

    /// Deform the body by updating the displacement x = X + u
    /// and update the internal variables with the new deformation and the
    /// time step increment
    void update_internal_variables(vector const& u, double const time_step_size = 0.0);

    /// Update the internal variables if converged, otherwise revert back
    /// for next attempted load increment
    void save_internal_variables(bool const have_converged);

    /// Constant access to the sub-meshes
    [[nodiscard]] std::vector<fem_submesh> const& meshes() const noexcept { return submeshes; }

    /// Non-const access to the sub-meshes
    [[nodiscard]] auto& meshes() noexcept { return submeshes; }

    [[nodiscard]] auto const& dirichlet_boundaries() const { return displacement_bcs; }

    [[nodiscard]] auto const& nonfollower_boundaries() const { return nonfollower_loads; }

    /// Gathers the time history for each boundary condition and
    /// returns a sorted vector which may contain traces of duplicates.
    /// \sa adaptive_time_step
    [[nodiscard]] std::vector<double> time_history() const;

    /// Provide const access to the discretised geometry for this mesh
    [[nodiscard]] auto const& geometry() const { return *coordinates; }

    /// Write out results to file
    void write(std::int64_t const time_step, double const current_time);

protected:
    void check_boundary_conditions(json const& boundary_data) const;

    void allocate_boundary_conditions(json const& boundary_data, basic_mesh const& reference_mesh);

    void allocate_displacement_boundary(json const& boundary, basic_mesh const& reference_mesh);

    [[nodiscard]] bool is_nonfollower_load(std::string const& boundary_type) const;

protected:
    std::shared_ptr<material_coordinates> coordinates;

    /// Meshes that contain individual element types
    std::vector<fem_submesh> submeshes;

    /// Displacement boundary conditions
    std::map<std::string, std::vector<dirichlet>> displacement_bcs;

    /// Nonfollower (force) boundary conditions
    std::map<std::string, nonfollower_load_boundary> nonfollower_loads;

    /// Nodal reaction forces
    vector reaction_forces;

    std::unordered_map<std::string, int> const dof_table = {{"x", 0}, {"y", 1}, {"z", 2}};

    /// This time step is taken from
    /// "Time[Period][Increments][Initial]" in the
    /// input file.  It is used in the boundary class
    /// to generate cyclic loading for example. This
    /// ensures the compatibility between user
    /// defined and sinusoidal boundary conditions.
    double generate_time_step;

    /// File output handle
    // std::unique_ptr<io::file_output> writer;
};
}
