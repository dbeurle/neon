
#pragma once

#include "mesh/boundary/neumann.hpp"

#include "newton_convection.hpp"

namespace neon
{
class basic_submesh;

namespace diffusion
{
/**
 * heat_flux is a group of the same elements on the same boundary. These
 * elements are responsible for computing their elemental right hand side contributions
 * with the corresponding shape function.  These are required to be stored in a parent
 * container with the other groups from the collective boundary \sa SurfaceBoundary
 */
using heat_flux = surface_load<surface_interpolation>;

// using heat_generation = volume_load<volume_interpolation>;

/* boundary_mesh contains the boundary conditions and meshes which contribute to
 * the external load vector.  This can include flux boundary conditions and Newton
 * convection type boundaries.  Each element group has an entry in the vector
 */
class boundary_mesh
{
public:
    explicit boundary_mesh(std::shared_ptr<material_coordinates>& material_coordinates,
                           std::vector<basic_submesh> const& submeshes,
                           json const& boundary,
                           json const& mesh_data);

    /** @return the boundaries which contribute only to the load vector */
    [[nodiscard]] auto const& load_interface() const { return load_boundaries; }

    /** @return the boundaries which contribute to the stiffness and the load vector */
    [[nodiscard]] auto const& stiffness_load_interface() const { return stiffness_load_boundaries; }

protected:
    std::vector<heat_flux> load_boundaries;

    std::vector<boundary::newton_convection> stiffness_load_boundaries;
};
}
}
