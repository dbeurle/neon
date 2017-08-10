
#pragma once

#include "interpolations/SurfaceInterpolation.hpp"
#include "mesh/SubMesh.hpp"
#include "mesh/solid/MaterialCoordinates.hpp"
#include "mesh/solid/boundary/Boundary.hpp"
#include "numeric/DenseTypes.hpp"

#include <memory>

namespace neon::solid
{
/**
 * NonfollowerLoad is a base class for the class of loads which do not follow
 * the body throughout the simulation.  Every non-follower load will contribute
 * to the external force vector, whether from a volume or a surface load.
 */
class NonFollowerLoad : public Boundary
{
public:
    explicit NonFollowerLoad(std::vector<List> const& nodal_connectivity,
                             std::shared_ptr<MaterialCoordinates>& material_coordinates,
                             double const prescribed_load,
                             bool const is_load_ramped,
                             int const dof_offset,
                             int const nodal_dofs = 3);

    /** @return an element external force vector for a given element */
    virtual std::tuple<List const&, Vector> external_force(
        int const element, double const load_factor) const = 0;

    auto elements() const { return nodal_connectivity.size(); }

protected:
    std::vector<List> nodal_connectivity;
    std::vector<List> dof_list;

    std::shared_ptr<MaterialCoordinates> material_coordinates;
};

/**
 * Traction is a non-follower load that has a surface interpolation and
 * computes the element external load vector contribution to the system of
 * equations \sa NonFollowerLoad
 */
class Traction : public NonFollowerLoad
{
public:
    /**
     * Construct the object by forwarding the constructor arguments for the
     * parent via a variadic template.  The correct constructor arguments
     * are given by \sa NonFollowerLoad
     * @param sf Unique pointer to a surface interpolation
     * @param args See \sa NonFollowerLoad constructor arguments
     */
    template <typename... NonFollowerLoadArgs>
    explicit Traction(std::unique_ptr<SurfaceInterpolation>&& sf,
                      NonFollowerLoadArgs... args)
        : NonFollowerLoad(args...), sf(std::move(sf))
    {
    }

    virtual std::tuple<List const&, Vector> external_force(
        int const element, double const load_factor) const override;

protected:
    std::unique_ptr<SurfaceInterpolation> sf;
};

/**
 * NonFollowerLoadBoundary contains the boundary conditions which contributed to
 * the external force vector.  This includes tractions, pressures and nodal
 * forces defined in the initial configuration
 */
class NonFollowerLoadBoundary
{
public:
    explicit NonFollowerLoadBoundary(std::shared_ptr<MaterialCoordinates>& material_coordinates,
                                     std::vector<SubMesh> const& submeshes,
                                     int const dof_offset,
                                     double const prescribed_load,
                                     Json::Value const& simulation_data);

    auto const& boundaries() const { return nf_loads; }

protected:
    std::vector<Traction> nf_loads;
};
}
