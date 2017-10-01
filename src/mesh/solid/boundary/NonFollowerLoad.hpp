
#pragma once

#include "mesh/common/Neumann.hpp"

#include "interpolations/InterpolationFactory.hpp"
#include "mesh/MaterialCoordinates.hpp"
#include "mesh/SubMesh.hpp"
#include "numeric/DenseTypes.hpp"

namespace neon::solid
{
/**
 * Traction is a non-follower load that has a surface interpolation and
 * computes the element external load vector contribution to the system of
 * equations \sa NonFollowerLoad
 */
using Traction = SurfaceLoad<SurfaceInterpolation>;

/**
 * Pressure is responsible for computing the resulting force addition from a
 * uniform pressure load that does not follow the deformation of the body.
 */
class Pressure : public Traction
{
public:
    Pressure(std::unique_ptr<SurfaceInterpolation>&& sf,
             std::vector<List> const& nodal_connectivity,
             std::shared_ptr<MaterialCoordinates>& material_coordinates,
             Json::Value const& time_history,
             Json::Value const& load_history,
             int const nodal_dofs);

    virtual std::tuple<List const&, Vector> external_force(int const element,
                                                           double const load_factor) const override;
};

/**
 * NonFollowerLoadBoundary contains the boundary conditions which contribute to
 * the external force vector.  This can include tractions, pressures, nodal
 * forces and volume forces computed in the initial configuration
 */
class NonFollowerLoadBoundary
{
public:
    explicit NonFollowerLoadBoundary(std::shared_ptr<MaterialCoordinates>& material_coordinates,
                                     std::vector<SubMesh> const& submeshes,
                                     Json::Value const& times,
                                     Json::Value const& loads,
                                     Json::Value const& simulation_data,
                                     int const dof_offset);

    auto const& boundaries() const { return surface_loads; }

protected:
    std::vector<Traction> surface_loads;
};
}
