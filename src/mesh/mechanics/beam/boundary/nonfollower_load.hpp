
#pragma once

/// @file

#include "interpolations/shape_function.hpp"
#include "mesh/basic_submesh.hpp"
#include "mesh/material_coordinates.hpp"

#include "mesh/boundary/nodal_value.hpp"
#include "couple.hpp"
#include "distributed_load.hpp"

#include <variant>

namespace neon::mechanics::beam
{
/// nonfollower_load_boundary contains the boundary conditions which contribute to
/// the external force vector.  This can include couples, distributed loads, nodal
/// forces computed in the initial configuration
class nonfollower_load_boundary
{
public:
    /// Specifying allowable nonfollower loads types
    using value_types = std::variant<couple, distributed_load>;

public:
    explicit nonfollower_load_boundary(std::shared_ptr<material_coordinates>& coordinates,
                                       std::vector<basic_submesh> const& submeshes,
                                       json const& simulation_data,
                                       json const& boundary_data,
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
