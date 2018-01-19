
#pragma once

#include "mesh/generic/Neumann.hpp"

#include "interpolations/InterpolationFactory.hpp"
#include "mesh/Submesh.hpp"

#include <variant>

namespace neon::mechanical::plane
{
/**
 * traction is a non-follower load that has a surface interpolation and
 * computes the element external load vector contribution to the system of
 * equations
 * \sa NonFollowerLoad
 */
using traction = SurfaceLoad<LineInterpolation>;

/**
 * body_force is a non-follower load that has a volume interpolation and
 * computes the element external load vector contribution to the system of
 * equations
 * \sa NonFollowerLoad
 */
using body_force = VolumeLoad<SurfaceInterpolation>;

/**
 * NonFollowerLoadBoundary contains the boundary conditions which contribute to
 * the external force vector.  This can include tractions, pressures, nodal
 * forces and volume forces computed in the initial configuration
 *
 * \sa traction
 * \sa Pressure
 * \sa body_force
 */
class NonFollowerLoadBoundary
{
public:
    /** Specifying the allowable nonfollower loads */
    using boundary_type = std::vector<std::variant<traction, body_force>>;

public:
    explicit NonFollowerLoadBoundary(std::shared_ptr<MaterialCoordinates>& material_coordinates,
                                     std::vector<Submesh> const& submeshes,
                                     Json::Value const& simulation_data,
                                     Json::Value const& boundary,
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
    std::array<std::pair<bool, boundary_type>, 2> nonfollower_load;
};
}
