
#pragma once

#include "mesh/boundary/dirichlet.hpp"
#include "mesh/mechanics/beam/boundary/nonfollower_load.hpp"
#include "mesh/mechanics/beam/submesh.hpp"
#include "io/file_output.hpp"

#include <map>

namespace neon
{
class basic_mesh;

namespace mechanics::beam
{
/// mesh is the discretised space for C0 beam finite elements.  It handles
/// the submeshes for each structural section (profile, material etc) and
/// is responsible for populating the submeshes and boundary conditions.
class mesh
{
public:
    using internal_variable_type = submesh::internal_variable_type;

    /// Alias traits to submesh
    using traits = submesh::traits;

public:
    mesh(basic_mesh const& basic_mesh,
         json const& material_data,
         json const& simulation_data,
         double const generate_time_step,
         std::map<std::string, std::unique_ptr<geometry::profile>> const& profile_store);

    /// The number of active degrees of freedom in this mesh
    [[nodiscard]] auto active_dofs() const noexcept
    {
        return traits::dofs_per_node * coordinates->size();
    }

    /// Checks the boundary conditions and constitutive model to ensure
    /// resulting matrix from this mesh is symmetric.
    [[nodiscard]] bool is_symmetric() const;

    /// Deform the body by updating the displacement x = X + u
    /// and update the internal variables with the new deformation and the
    /// time step increment
    void update_internal_variables(vector const& displacement_rotation,
                                   double const time_step_size = 0.0);

    /// Constant access to the sub-meshes
    [[nodiscard]] auto const& meshes() const noexcept { return submeshes; }

    /// Non-const access to the sub-meshes
    [[nodiscard]] auto& meshes() noexcept { return submeshes; }

    [[nodiscard]] auto const& dirichlet_boundaries() const noexcept
    {
        return m_dirichlet_boundaries;
    }

    /// Boundary conditions associated with consistent (element) loads or
    /// nodal based loads
    /// \sa nonfollower_load_boundary
    /// \return vector of nonfollower boundaries
    [[nodiscard]] auto const& nonfollower_boundaries() const noexcept { return nonfollower_loads; }

    /// Gathers the time history for each boundary condition and
    /// returns a sorted vector which may contain traces of duplicates.
    /// \sa adaptive_time_step
    [[nodiscard]] std::vector<double> time_history() const;

    /// Provide const access to the discretised geometry for this mesh
    [[nodiscard]] auto const& geometry() const { return *coordinates; }

    /// Write out results to file
    void write(std::int32_t const time_step, double const current_time);

protected:
    void check_boundary_conditions(json const& boundary_data) const;

    void allocate_variable_names();

    void allocate_boundary_conditions(json const& boundary_data, basic_mesh const& basic_mesh);

    void allocate_dirichlet_boundary(std::string const& boundary_type,
                                     json const& boundary,
                                     basic_mesh const& basic_mesh);

    [[nodiscard]] bool is_nonfollower_load(std::string const& boundary_type) const;

protected:
    std::shared_ptr<material_coordinates> coordinates;

    std::vector<submesh> submeshes;

    /// Displacement boundary conditions
    std::map<std::string, std::vector<dirichlet>> m_dirichlet_boundaries;

    /// Nonfollower (force and moment) boundary conditions
    std::map<std::string, nonfollower_load_boundary> nonfollower_loads;

    std::unordered_map<std::string, int> const dof_table{{"x", 0}, {"y", 1}, {"z", 2}};

    /// Profiles available for the submeshes
    std::map<std::string, geometry::profile> profiles;

    /// Displacement and rotation vector (u1, u2, u3, r1, r2, r3)
    vector displacement, rotation;

    /// This time step is taken from
    /// "Time[Period][Increments][Initial]" in the
    /// input file.  It is used in the boundary class
    /// to generate cyclic loading for example. This
    /// ensures the compatibility between user
    /// defined and sinusoidal boundary conditions.
    double generate_time_step;

    /// File output handle
    std::unique_ptr<io::file_output> writer;

    /// Output variables
    std::vector<variable::types> output_variables;
};
}
}
