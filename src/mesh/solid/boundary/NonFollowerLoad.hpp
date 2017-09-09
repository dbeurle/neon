
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
                                     int const dof_offset)
    {
        for (auto const& mesh : submeshes)
        {
            surface_loads.emplace_back(make_surface_interpolation(mesh.topology(), simulation_data),
                                       mesh.connectivities(),
                                       material_coordinates,
                                       times,
                                       loads,
                                       dof_offset,
                                       3);
        }
    }

    auto const& boundaries() const { return surface_loads; }

protected:
    std::vector<Traction> surface_loads;
};
}
