
#pragma once

#include "body_force.hpp"
#include "pressure.hpp"
#include "traction.hpp"

#include "interpolations/shape_function.hpp"
#include "mesh/basic_submesh.hpp"

#include <variant>

namespace neon::mechanical::solid
{
/**
 * NonFollowerLoadBoundary contains the boundary conditions which contribute to
 * the external force vector.  This can include tractions, pressures, nodal
 * forces and volume forces computed in the initial configuration
 *
 * \sa traction
 * \sa pressure
 * \sa body_force
 */
class NonFollowerLoadBoundary
{
public:
    using value_types = std::variant<traction, pressure, body_force>;

    /** Specifying the allowable nonfollower loads */
    using BoundaryMeshes = std::vector<std::variant<traction, pressure, body_force>>;

public:
    explicit NonFollowerLoadBoundary(std::shared_ptr<material_coordinates>& material_coordinates,
                                     std::vector<basic_submesh> const& submeshes,
                                     json const& simulation_data,
                                     json const& boundary,
                                     std::unordered_map<std::string, int> const& dof_table);

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
    std::array<std::pair<bool, BoundaryMeshes>, 3> nonfollower_load;
};
}
