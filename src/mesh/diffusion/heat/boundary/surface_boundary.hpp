
#pragma once

#include "mesh/boundary/neumann.hpp"

#include "newton_convection.hpp"

namespace neon
{
class basic_submesh;

namespace diffusion
{
// using heat_generation = volume_load<volume_interpolation>;

/// boundary_mesh contains the boundary conditions and meshes which contribute to
/// the external load vector.  This can include flux boundary conditions and Newton
/// convection type boundaries.  Each element group has an entry in the vector
class boundary_mesh
{
public:
    explicit boundary_mesh(std::shared_ptr<material_coordinates>& material_coordinates,
                           std::vector<basic_submesh> const& submeshes,
                           json const& boundary_data,
                           json const& mesh_data);

    /// \return the boundaries which contribute only to the load vector
    [[nodiscard]] auto const& load_interface() const { return load_boundaries; }

    /// \return the boundaries which contribute to the stiffness and the load vector
    [[nodiscard]] auto const& stiffness_load_interface() const { return stiffness_load_boundaries; }

protected:
    std::vector<boundary::heat_flux> load_boundaries;

    std::vector<boundary::newton_convection> stiffness_load_boundaries;
};
}
}
