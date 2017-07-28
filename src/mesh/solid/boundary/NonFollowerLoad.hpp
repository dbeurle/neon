
#pragma once

#include "interpolations/SurfaceInterpolation.hpp"
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

class Traction : public NonFollowerLoad
{
public:
    Traction(std::vector<List> const& nodal_connectivity,
             std::shared_ptr<MaterialCoordinates>&& material_coordinates,
             double const prescribed_load,
             bool const is_load_ramped,
             int const dof_offset,
             std::unique_ptr<SurfaceInterpolation>&& sf);

    virtual std::tuple<List const&, Vector> external_force(
        int const element, double const load_factor) const override;

protected:
    std::unique_ptr<SurfaceInterpolation> sf;
};
}
