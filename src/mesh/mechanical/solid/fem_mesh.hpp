
#pragma once

#include "mesh/generic/Dirichlet.hpp"
#include "mesh/mechanical/solid/boundary/NonFollowerLoad.hpp"
#include "mesh/mechanical/solid/fem_submesh.hpp"

#include <map>

namespace neon
{
class basic_mesh;

namespace mechanical::solid
{
class fem_mesh
{
public:
    using internal_variable_type = fem_submesh::internal_variable_type;

public:
    fem_mesh(basic_mesh const& basic_mesh,
             json const& material_data,
             json const& simulation_data,
             double const generate_time_step);

    /** The number of active degrees of freedom in this mesh */
    [[nodiscard]] auto active_dofs() const { return 3 * mesh_coordinates->size(); }

    /**
     * Checks the boundary conditions and constitutive model to ensure
     * resulting matrix from this mesh is symmetric.  \sa LinearSolver
     */
    [[nodiscard]] bool is_symmetric() const;

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
    [[nodiscard]] std::vector<fem_submesh> const& meshes() const { return submeshes; }

        /** Non-const access to the sub-meshes */
        [[nodiscard]] std::vector<fem_submesh>& meshes()
    {
        return submeshes;
    }

    [[nodiscard]] auto const& displacement_boundaries() const { return displacement_bcs; }

    [[nodiscard]] auto const& nonfollower_load_boundaries() const { return nonfollower_loads; }

    /**
     * Gathers the time history for each boundary condition and
     * returns a sorted vector which may contain traces of duplicates.
     *
     * \sa adaptive_time_step
     */
    [[nodiscard]] std::vector<double> time_history() const;

    [[deprecated]][[nodiscard]] auto const& coordinates() const { return *mesh_coordinates; }

    /** Provide const access to the discretised geometry for this mesh */
    [[nodiscard]] auto const& geometry() const { return *mesh_coordinates; }

protected:
    void check_boundary_conditions(json const& boundary_data) const;

    void allocate_boundary_conditions(json const& boundary_data, basic_mesh const& basic_mesh);

    void allocate_displacement_boundary(json const& boundary, basic_mesh const& basic_mesh);

    [[nodiscard]] bool is_nonfollower_load(std::string const& boundary_type) const;

    /** Collapse the nodal connectivity arrays from the submesh for a node list */
    [[nodiscard]] local_indices filter_dof_list(std::vector<basic_submesh> const& boundary_mesh) const;

protected:
    /**
     * This time step is taken from "Time[Period][Increments][Initial]" in the input file.
     * It is used in the boundary class to generate cyclic loading for example. This ensures the
     compatibility between user defined and sinusoidal boundary conditions.
     */
    double generate_time_step;

    std::shared_ptr<material_coordinates> mesh_coordinates;

    std::vector<fem_submesh> submeshes;

    // Boundary conditions for this mesh
    std::map<std::string, std::vector<Dirichlet>> displacement_bcs;
    std::map<std::string, NonFollowerLoadBoundary> nonfollower_loads;

    std::unordered_map<std::string, int> const dof_table = {{"x", 0}, {"y", 1}, {"z", 2}};
};
}
}
