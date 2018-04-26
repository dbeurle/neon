
#pragma once

#include "mesh/generic/neumann.hpp"

#include "interpolations/interpolation_factory.hpp"
#include "mesh/basic_submesh.hpp"

#include <variant>

namespace neon::mechanical::plane
{
/// traction is a non-follower load that has a surface interpolation and
/// computes the element external load vector contribution to the system of
/// equations
/// \sa NonFollowerLoad
using traction = surface_load<line_interpolation>;

/// body_force is a non-follower load that has a volume interpolation and
/// computes the element external load vector contribution to the system of
/// equations
/// \sa NonFollowerLoad
using body_force = volume_load<surface_interpolation>;

/// nonfollower_load_boundary contains the boundary conditions which contribute to
/// the external force vector.  This can include tractions, pressures, nodal
/// forces and volume forces computed in the initial configuration
///
/// \sa traction
/// \sa Pressure
/// \sa body_force
class nonfollower_load_boundary
{
public:
<<<<<<< HEAD
    /// Specify the allowable nonfollower loads
=======
    /// Specifying allowable nonfollower loads types
>>>>>>> 3827a3fbd98c45ebe38532028944988f58b9737d
    using boundary_type = std::vector<std::variant<traction, body_force>>;

public:
    explicit nonfollower_load_boundary(std::shared_ptr<material_coordinates>& material_coordinates,
                                       std::vector<basic_submesh> const& submeshes,
                                       json const& simulation_data,
                                       json const& boundary,
                                       std::unordered_map<std::string, int> const& dof_table,
                                       double const generate_time_step);

    /**
     * Provides access to an array with three elements.  Each element is a pair,
     * with the first corresponds to an active DoF flag and the second is a
     * variant type of the allowable boundary conditions.  This can be accessed
     * using

           for (auto const& [name, nonfollower_load] : mesh.nonfollower_loads())
           {
               for (auto const& [is_dof_active, boundary_condition] : nonfollower_load.interface())
               {
                   if (!is_dof_active) continue;

                   std::visit([](auto const& mesh) {
                       for (auto element = 0; element < mesh.elements(); ++element)
                       {
                            External force assembly
                       }
                   }, boundary_condition);
               }
           }

     */
    auto const& interface() const { return nonfollower_load; }

protected:
    std::array<std::pair<bool, boundary_type>, 2> nonfollower_load;
};
}
