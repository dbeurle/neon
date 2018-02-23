
#pragma once

#include "vector_contribution.hpp"

#include "geometry/Projection.hpp"
#include "mesh/material_coordinates.hpp"

#include <memory>

namespace neon
{
/**
 * neumann is a base class for Neumann (derivative) type boundary conditions.
 * This includes the nodal connectivities and degrees of freedom lists.  Derived
 * classes must implement shape functions and the appropriate finite element
 * approximation for the given problem
 */
class neumann : public vector_contribution
{
public:
    explicit neumann(indices nodal_connectivity,
                     indices dof_list,
                     std::shared_ptr<material_coordinates>& material_coordinates,
                     json const& times,
                     json const& loads);

    explicit neumann(indices nodal_connectivity,
                     indices dof_list,
                     std::shared_ptr<material_coordinates>& material_coordinates,
                     json const& boundary,
                     std::string const& name,
                     double const generate_time_step);

    [[nodiscard]] auto elements() const { return nodal_connectivity.cols(); }

protected:
    indices nodal_connectivity, dof_list;

    std::shared_ptr<material_coordinates> coordinates;
};

template <typename SurfaceInterpolation_Tp>
class surface_load : public neumann
{
public:
    explicit surface_load(std::unique_ptr<SurfaceInterpolation_Tp>&& sf,
                          indices nodal_connectivity,
                          indices dof_list,
                          std::shared_ptr<material_coordinates>& coordinates,
                          json const& time_history,
                          json const& load_history)
        : neumann(nodal_connectivity, dof_list, coordinates, time_history, load_history),
          sf(std::move(sf))
    {
    }

    explicit surface_load(std::unique_ptr<SurfaceInterpolation_Tp>&& sf,
                          indices nodal_connectivity,
                          indices dof_list,
                          std::shared_ptr<material_coordinates>& coordinates,
                          json const& boundary,
                          std::string const& name,
                          double const generate_time_step)
        : neumann(nodal_connectivity, dof_list, coordinates, boundary, name, generate_time_step),
          sf(std::move(sf))
    {
    }

    /**
     * Compute the external force due to a neumann type boundary condition.
     * This computes the following integral on a boundary element
     */
    virtual std::pair<index_view, vector> external_force(std::int64_t const element,
                                                         double const load_factor) const override
    {
        auto const X = geometry::project_to_plane(coordinates->initial_configuration(
            nodal_connectivity(Eigen::placeholders::all, element)));

        // Perform the computation of the external load vector
        auto const f_ext = sf->quadrature().integrate(vector::Zero(X.cols()).eval(),
                                                      [&](auto const& femval, auto const& l) -> vector {
                                                          auto const& [N, dN] = femval;

                                                          auto const j = (X * dN).determinant();

                                                          return N * j;
                                                      });

        return {dof_list(Eigen::placeholders::all, element),
                interpolate_prescribed_load(load_factor) * f_ext};
    }

protected:
    std::unique_ptr<SurfaceInterpolation_Tp> sf;
};

template <typename VolumeInterpolation_Tp>
class volume_load : public neumann
{
public:
    explicit volume_load(std::unique_ptr<VolumeInterpolation_Tp>&& sf,
                         indices nodal_connectivity,
                         indices dof_list,
                         std::shared_ptr<material_coordinates>& coordinates,
                         json const& time_history,
                         json const& load_history)
        : neumann(nodal_connectivity, dof_list, coordinates, time_history, load_history),
          sf(std::move(sf))
    {
    }

    explicit volume_load(std::unique_ptr<VolumeInterpolation_Tp>&& sf,
                         indices nodal_connectivity,
                         indices dof_list,
                         std::shared_ptr<material_coordinates>& coordinates,
                         json const& boundary,
                         std::string const& name,
                         double const generate_time_step)
        : neumann(nodal_connectivity, dof_list, coordinates, boundary, name, generate_time_step),
          sf(std::move(sf))
    {
    }

    std::pair<index_view, vector> external_force(std::int64_t const element,
                                                 double const load_factor) const override
    {
        auto const X = coordinates->initial_configuration(
            nodal_connectivity(Eigen::placeholders::all, element));

        // Perform the computation of the external load vector
        auto const f_ext = sf->quadrature().integrate(vector::Zero(X.cols()).eval(),
                                                      [&](auto const& femval, auto const& l) -> vector {
                                                          auto const& [N, dN] = femval;

                                                          auto const j = (X * dN).determinant();

                                                          return N * j;
                                                      });

        return {dof_list(Eigen::placeholders::all, element),
                interpolate_prescribed_load(load_factor) * f_ext};
    }

protected:
    std::unique_ptr<VolumeInterpolation_Tp> sf;
};
}
