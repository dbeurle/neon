
#pragma once

/// @file

#include "mesh/boundary/nodal_value.hpp"
#include "mesh/mechanics/plane/boundary/traction.hpp"
#include "mesh/mechanics/plane/boundary/body_force.hpp"

#include "mesh/basic_submesh.hpp"

#include <variant>

namespace neon::mechanics::plane
{
/// nonfollower_load_boundary contains the boundary conditions which contribute to
/// the external force vector.  This can include tractions, pressures, nodal
/// forces and volume forces computed in the initial configuration
/// \sa traction
/// \sa body_force
class nonfollower_load_boundary
{
public:
    /// Specifying allowable nonfollower loads types
    using value_types = std::variant<traction, body_force>;

public:
    explicit nonfollower_load_boundary(std::shared_ptr<material_coordinates>& coordinates,
                                       std::vector<basic_submesh> const& submeshes,
                                       json const& simulation_data,
                                       json const& boundary,
                                       std::unordered_map<std::string, int> const& dof_table,
                                       double const generate_time_step);

    /// Provides const access to a vector of variants containing the types
    /// that represent natural boundary conditions on the surface of the object
    auto const& natural_interface() const noexcept { return boundary_meshes; }

    auto const& nodal_interface() const noexcept { return nodal_values; }

protected:
    std::vector<value_types> boundary_meshes;

    std::vector<nodal_value> nodal_values;
};
}
