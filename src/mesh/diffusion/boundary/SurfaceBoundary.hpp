
#pragma once

#include "mesh/generic/Neumann.hpp"

#include "NewtonConvection.hpp"

#include "interpolations/InterpolationFactory.hpp"
#include "mesh/diffusion/Submesh.hpp"

namespace neon::diffusion
{
/**
 * SurfaceFlux is a group of the same elements on the same boundary. These
 * elements are responsible for computing their elemental right hand side contributions
 * with the corresponding shape function.  These are required to be stored in a parent
 * container with the other groups from the collective boundary \sa SurfaceBoundary
 */
using SurfaceFlux = SurfaceLoad<SurfaceInterpolation>;

/* SurfaceBoundary contains the boundary conditions which contribute to
 * the external load vector.  This can include flux boundary conditions and Newton
 * convection type boundaries.  Each element group has an entry in the vector
 */
class SurfaceBoundary
{
public:
    explicit SurfaceBoundary(std::shared_ptr<MaterialCoordinates>& material_coordinates,
                             std::vector<Submesh> const& submeshes,
                             Json::Value const& times,
                             Json::Value const& loads,
                             Json::Value const& simulation_data)
    {
        for (auto const& mesh : submeshes)
        {
            surface_forces.emplace_back(make_surface_interpolation(mesh.topology(), simulation_data),
                                        mesh.connectivities(),
                                        mesh.connectivities(),
                                        material_coordinates,
                                        times,
                                        loads);
        }
    }

    auto const& load_interface() const { return surface_forces; }

    auto const& stiffness_interface() const { return surface_stiffness_loads; }

protected:
    std::vector<SurfaceFlux> surface_forces;

    std::vector<NewtonConvection> surface_stiffness_loads;
};
}
