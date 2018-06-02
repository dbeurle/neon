
#pragma once

#include "interpolations/shape_function.hpp"
#include "mesh/basic_submesh.hpp"

#include "mesh/generic/nodal_value.hpp"
#include "body_force.hpp"
#include "pressure.hpp"
#include "traction.hpp"

#include <variant>

namespace neon::mechanical::solid
{
/// nonfollower_load_boundary contains the boundary conditions which contribute to
/// the external force vector.  This can include tractions, pressures, nodal
/// forces and volume forces computed in the initial configuration
/// \sa traction
/// \sa pressure
/// \sa body_force
class nonfollower_load_boundary
{
public:
    /// Specifying allowable nonfollower loads types
    using value_types = std::variant<traction, pressure, body_force>;

public:
    explicit nonfollower_load_boundary(std::shared_ptr<material_coordinates>& material_coordinates,
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
