
#pragma once

#include "mesh/generic/Dirichlet.hpp"
#include "mesh/mechanical/plane/boundary/NonFollowerLoad.hpp"
#include "mesh/mechanical/plane/femSubmesh.hpp"

#include <map>
#include <unordered_map>

namespace neon
{
class BasicMesh;

namespace mechanical::plane
{
class femMesh
{
public:
    using internal_variable_type = InternalVariables;

public:
    femMesh(BasicMesh const& basic_mesh,
            Json::Value const& material_data,
            Json::Value const& simulation_data);

    /** The number of active degrees of freedom in this mesh */
    [[nodiscard]] auto active_dofs() const { return 2 * material_coordinates->size(); }

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
    [[nodiscard]] std::vector<femSubmesh> const& meshes() const { return submeshes; }

        /** Mutable access to the sub-meshes */
        [[nodiscard]] std::vector<femSubmesh>& meshes()
    {
        return submeshes;
    }

    [[nodiscard]] auto const& displacement_boundaries() const { return displacement_bcs; }

    [[nodiscard]] auto const& nonfollower_load_boundaries() const { return nonfollower_loads; }

    /**
     * Gathers the time history for each boundary condition and
     * returns a sorted vector which may contain traces of duplicates.
     *
     * \sa AdaptiveLoadStep
     */
    [[nodiscard]] std::vector<double> time_history() const;

    [[nodiscard]] auto const& coordinates() const { return *(material_coordinates.get()); }

protected:
    void check_boundary_conditions(Json::Value const& boundary_data) const;

    void allocate_boundary_conditions(Json::Value const& boundary_data, BasicMesh const& basic_mesh);

    void allocate_displacement_boundary(Json::Value const& boundary, BasicMesh const& basic_mesh);

    [[nodiscard]] bool is_nonfollower_load(std::string const& boundary_type) const;

    /** Collapse the nodal connectivity arrays from the submesh for a node list */
    [[nodiscard]] List filter_dof_list(std::vector<Submesh> const& boundary_mesh) const;

protected:
    std::shared_ptr<MaterialCoordinates> material_coordinates;

    std::vector<femSubmesh> submeshes;

    // Boundary conditions for this mesh
    std::map<std::string, std::vector<Dirichlet>> displacement_bcs;
    std::map<std::string, NonFollowerLoadBoundary> nonfollower_loads;

    std::unordered_map<std::string, int> const dof_table = {{"x", 0}, {"y", 1}};
};
}
}
