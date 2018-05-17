
#pragma once

#include "mesh/generic/dirichlet.hpp"
#include "mesh/mechanical/plane/boundary/nonfollower_load.hpp"
#include "mesh/mechanical/plane/fem_submesh.hpp"

#include <map>
#include <vector>
#include <unordered_map>

namespace neon
{
class basic_mesh;

namespace mechanical::plane
{
class fem_mesh
{
public:
    using internal_variable_type = internal_variables_t;

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

    void update_internal_forces(vector const& fint) { internal_forces = fint; }

    vector nodal_reaction_forces() const { return -internal_forces; }

    /// Deform the body by updating the displacement x = X + u
    /// and update the internal variables with the new deformation and the
    /// time step increment
    void update_internal_variables(vector const& u, double const time_step_size = 0.0);

    /// Update the internal variables if converged, otherwise revert back
    /// for next attempted load increment
    void save_internal_variables(bool const have_converged);

    /// Constant access to the sub-meshes
    [[nodiscard]] std::vector<fem_submesh> const& meshes() const noexcept { return submeshes; }

    /// Mutable access to the sub-meshes
    [[nodiscard]] std::vector<fem_submesh>& meshes() noexcept { return submeshes; }

    [[nodiscard]] auto const& dirichlet_boundaries() const { return displacement_bcs; }

    [[nodiscard]] auto const& nonfollower_load_boundaries() const { return nonfollower_loads; }

    /// Gathers the time history for each boundary condition and
    /// returns a sorted vector which may contain duplicated entries.
    /// \sa adaptive_time_step
    [[nodiscard]] std::vector<double> time_history() const;

    [[nodiscard]] auto const& geometry() const { return *coordinates; }

protected:
    void check_boundary_conditions(json const& boundary_data) const;

    void allocate_boundary_conditions(json const& boundary_data, basic_mesh const& basic_mesh);

    void allocate_displacement_boundary(json const& boundary, basic_mesh const& basic_mesh);

    [[nodiscard]] bool is_nonfollower_load(std::string const& boundary_type) const;

protected:
    std::shared_ptr<material_coordinates> coordinates;

    std::vector<fem_submesh> submeshes;

    /// Displacement boundaries
    std::map<std::string, std::vector<dirichlet>> displacement_bcs;
    /// Nonfollower boundaries
    std::map<std::string, nonfollower_load_boundary> nonfollower_loads;

    /// Internal nodal forces for reaction forces
    vector internal_forces;

    std::unordered_map<std::string, int> const dof_table = {{"x", 0}, {"y", 1}};

    /// This time step is taken from
    /// "Time[Period][Increments][Initial]" in the
    /// input file.  It is used in the boundary class
    /// to generate cyclic loading for example. This
    /// ensures the compatibility between user
    /// defined and sinusoidal boundary conditions.
    double generate_time_step;
};
}
}
