
#pragma once

#include "mesh/common/Neumann.hpp"

#include "interpolations/InterpolationFactory.hpp"
#include "mesh/SubMesh.hpp"

#include <variant>

namespace neon::solid
{
/**
 * Traction is a non-follower load that has a surface interpolation and
 * computes the element external load vector contribution to the system of
 * equations
 * \sa NonFollowerLoad
 */
using Traction = SurfaceLoad<SurfaceInterpolation>;

/**
 * BodyForce is a non-follower load that has a volume interpolation and
 * computes the element external load vector contribution to the system of
 * equations
 * \sa NonFollowerLoad
 */
using BodyForce = VolumeLoad<VolumeInterpolation>;

class Pressure : public Traction
{
public:
    explicit Pressure(std::unique_ptr<SurfaceInterpolation>&& sf,
                      std::vector<List> const& nodal_connectivity,
                      std::shared_ptr<MaterialCoordinates>& material_coordinates,
                      Json::Value const& time_history,
                      Json::Value const& load_history,
                      int const nodal_dofs);

    std::tuple<List const&, Vector> external_force(int const element,
                                                   double const load_factor) const override;
};

/**
 * NonFollowerLoadBoundary contains the boundary conditions which contribute to
 * the external force vector.  This can include tractions, pressures, nodal
 * forces and volume forces computed in the initial configuration
 *
 * \sa Traction
 * \sa Pressure
 * \sa BodyForce
 */
class NonFollowerLoadBoundary
{
public:
    /** Specifying the allowable nonfollower loads */
    using BoundaryMeshes = std::vector<std::variant<Traction, Pressure, BodyForce>>;

public:
    explicit NonFollowerLoadBoundary(std::shared_ptr<MaterialCoordinates>& material_coordinates,
                                     std::vector<SubMesh> const& submeshes,
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
    std::array<std::pair<bool, BoundaryMeshes>, 3> nonfollower_load;
};
}
