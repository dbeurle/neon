
#pragma once

#include "vector_contribution.hpp"

#include "geometry/Projection.hpp"
#include "mesh/DofAllocator.hpp"
#include "mesh/MaterialCoordinates.hpp"

#include <memory>

namespace neon::boundary
{
/**
 * neumann is a base class for neumann type boundary conditions.  This includes
 * the nodal connectivities and degrees of freedom lists.  Derived classes must
 * implement shape functions and the appropriate finite element approximation
 * for the given problem
 */
class neumann : public vector_contribution
{
public:
    explicit neumann(std::vector<std::vector<int64>> const& nodal_connectivity,
                     std::vector<std::vector<int64>> const& dof_list,
                     std::shared_ptr<MaterialCoordinates>& material_coordinates,
                     Json::Value const& times,
                     Json::Value const& loads);

    [[nodiscard]] auto elements() const { return nodal_connectivity.size(); }

protected:
    std::vector<std::vector<int64>> nodal_connectivity, dof_list;

    std::shared_ptr<MaterialCoordinates> material_coordinates;
};

template <typename SurfaceInterpolation_Tp>
class surface_load : public neumann
{
public:
    explicit surface_load(std::unique_ptr<SurfaceInterpolation_Tp>&& sf,
                          std::vector<std::vector<int64>> const& nodal_connectivity,
                          std::vector<std::vector<int64>> const& dof_list,
                          std::shared_ptr<MaterialCoordinates>& material_coordinates,
                          Json::Value const& time_history,
                          Json::Value const& load_history)
        : neumann(nodal_connectivity, dof_list, material_coordinates, time_history, load_history),
          sf(std::move(sf))
    {
    }

    /**
     * Compute the external force due to a Neumann type boundary condition.
     * This computes the following integral on a boundary element
     */
    virtual std::pair<std::vector<int64> const&, vector> external_force(
        int64 const element, double const load_factor) const override
    {
        auto const X = geometry::project_to_plane(
            material_coordinates->initial_configuration(nodal_connectivity[element]));

        // Perform the computation of the external load vector
        auto const f_ext = sf->quadrature().integrate(vector::Zero(X.cols()).eval(),
                                                      [&](auto const& femval, auto const& l) -> vector {
                                                          auto const& [N, dN] = femval;

                                                          auto const j = (X * dN).determinant();

                                                          return N * j;
                                                      });
        return {dof_list[element], interpolate_prescribed_load(load_factor) * f_ext};
    }

protected:
    std::unique_ptr<SurfaceInterpolation_Tp> sf;
};

template <typename VolumeInterpolation_Tp>
class volume_load : public neumann
{
public:
    explicit volume_load(std::unique_ptr<VolumeInterpolation_Tp>&& sf,
                         std::vector<std::vector<int64>> const& nodal_connectivity,
                         std::vector<std::vector<int64>> const& dof_list,
                         std::shared_ptr<MaterialCoordinates>& material_coordinates,
                         Json::Value const& time_history,
                         Json::Value const& load_history)
        : neumann(nodal_connectivity, dof_list, material_coordinates, time_history, load_history),
          sf(std::move(sf))
    {
    }

    std::pair<std::vector<int64> const&, vector> external_force(int64 const element,
                                                                double const load_factor) const override
    {
        auto const X = material_coordinates->initial_configuration(nodal_connectivity[element]);

        // Perform the computation of the external load vector
        auto const f_ext = sf->quadrature().integrate(vector::Zero(X.cols()).eval(),
                                                      [&](auto const& femval, auto const& l) -> vector {
                                                          auto const& [N, dN] = femval;

                                                          auto const j = (X * dN).determinant();

                                                          return N * j;
                                                      });
        return {dof_list.at(element), interpolate_prescribed_load(load_factor) * f_ext};
    }

protected:
    std::unique_ptr<VolumeInterpolation_Tp> sf;
};
}
