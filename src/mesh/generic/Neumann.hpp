
#pragma once

#include "VectorContribution.hpp"

#include "geometry/Projection.hpp"
#include "mesh/DofAllocator.hpp"
#include "mesh/MaterialCoordinates.hpp"

#include <memory>

namespace neon
{
/**
 * Neumann is a base class for Neumann type boundary conditions.  This includes
 * the nodal connectivities and degrees of freedom lists.  Derived classes must
 * implement shape functions and the appropriate finite element approximation
 * for the given problem
 */
class Neumann : public VectorContribution
{
public:
    explicit Neumann(std::vector<List> const& nodal_connectivity,
                     std::vector<List> const& dof_list,
                     std::shared_ptr<MaterialCoordinates>& material_coordinates,
                     json const& times,
                     json const& loads);

    [[nodiscard]] auto elements() const { return nodal_connectivity.size(); }

protected:
    std::vector<List> nodal_connectivity, dof_list;

    std::shared_ptr<MaterialCoordinates> material_coordinates;
};

template <typename SurfaceInterpolation_Tp>
class SurfaceLoad : public Neumann
{
public:
    explicit SurfaceLoad(std::unique_ptr<SurfaceInterpolation_Tp>&& sf,
                         std::vector<List> const& nodal_connectivity,
                         std::vector<List> const& dof_list,
                         std::shared_ptr<MaterialCoordinates>& material_coordinates,
                         json const& time_history,
                         json const& load_history)
        : Neumann(nodal_connectivity, dof_list, material_coordinates, time_history, load_history),
          sf(std::move(sf))
    {
    }

    /**
     * Compute the external force due to a Neumann type boundary condition.
     * This computes the following integral on a boundary element
     */
    virtual std::tuple<List const&, vector> external_force(int const element,
                                                           double const load_factor) const override
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
class VolumeLoad : public Neumann
{
public:
    explicit VolumeLoad(std::unique_ptr<VolumeInterpolation_Tp>&& sf,
                        std::vector<List> const& nodal_connectivity,
                        std::vector<List> const& dof_list,
                        std::shared_ptr<MaterialCoordinates>& material_coordinates,
                        json const& time_history,
                        json const& load_history)
        : Neumann(nodal_connectivity, dof_list, material_coordinates, time_history, load_history),
          sf(std::move(sf))
    {
    }

    std::tuple<List const&, vector> external_force(int const element,
                                                   double const load_factor) const override
    {
        auto const X = material_coordinates->initial_configuration(nodal_connectivity[element]);

        auto const f = interpolate_prescribed_load(load_factor);

        // Perform the computation of the external load vector
        auto const f_ext = sf->quadrature().integrate(vector::Zero(X.cols()).eval(),
                                                      [&](auto const& femval, auto const& l) -> vector {
                                                          auto const& [N, dN] = femval;

                                                          auto const j = (X * dN).determinant();

                                                          return f * N * j;
                                                      });
        return {dof_list.at(element), interpolate_prescribed_load(load_factor) * f_ext};
    }

protected:
    std::unique_ptr<VolumeInterpolation_Tp> sf;
};
}
