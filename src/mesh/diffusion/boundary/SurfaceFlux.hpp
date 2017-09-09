
#pragma once

#include "mesh/common/Neumann.hpp"

#include "interpolations/InterpolationFactory.hpp"
#include "mesh/diffusion/Submesh.hpp"

namespace neon::diffusion
{
/**
 * SurfaceFlux is a group of the same elements on the same boundary. These
 * elements are responsible for computing their elemental right hand side contributions
 * with the corresponding shape function.  These are required to be stored in a parent
 * container with the other groups from the collective boundary \sa SurfaceFluxBoundary
 */
using SurfaceFlux = SurfaceLoad<SurfaceInterpolation>;

/* SurfaceFluxBoundary contains the boundary conditions which contribute to
 * the external force vector.  This can include tractions, pressures, nodal
 * forces and volume forces computed in the initial configuration
 */
class SurfaceFluxBoundary
{
public:
    explicit SurfaceFluxBoundary(std::shared_ptr<MaterialCoordinates>& material_coordinates,
                                 std::vector<SubMesh> const& submeshes,
                                 Json::Value const& times,
                                 Json::Value const& loads,
                                 Json::Value const& simulation_data)
    {
        for (auto const& mesh : submeshes)
        {
            surface_loads.emplace_back(make_surface_interpolation(mesh.topology(), simulation_data),
                                       mesh.connectivities(),
                                       material_coordinates,
                                       times,
                                       loads,
                                       0,
                                       1);
        }
    }

    auto const& boundaries() const { return surface_loads; }

protected:
    std::vector<SurfaceFlux> surface_loads;
};
}
