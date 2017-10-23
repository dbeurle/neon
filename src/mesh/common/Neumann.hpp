
#pragma once

#include "Boundary.hpp"

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
class Neumann : public Boundary
{
public:
    explicit Neumann(std::vector<List> const& nodal_connectivity,
                     std::vector<List> const& dof_list,
                     std::shared_ptr<MaterialCoordinates>& material_coordinates,
                     Json::Value const& times,
                     Json::Value const& loads);

    auto elements() const { return nodal_connectivity.size(); }

    /** @return an element external force vector for a given element */
    virtual std::tuple<List const&, Vector> external_force(int const element,
                                                           double const load_factor) const = 0;

protected:
    std::vector<List> nodal_connectivity;
    std::vector<List> dof_list;

    std::shared_ptr<MaterialCoordinates> material_coordinates;
};

template <typename SurfaceInterpolation_Tp>
class SurfaceLoad : public Neumann
{
public:
    explicit SurfaceLoad(std::unique_ptr<SurfaceInterpolation_Tp>&& sf,
                         std::vector<List> const& nodal_connectivity,
                         std::shared_ptr<MaterialCoordinates>& material_coordinates,
                         Json::Value const& time_history,
                         Json::Value const& load_history,
                         int const dof_offset,
                         int const nodal_dofs)
        : Neumann(nodal_connectivity,
                  filter_dof_list(nodal_dofs, dof_offset, nodal_connectivity),
                  material_coordinates,
                  time_history,
                  load_history),
          sf(std::move(sf))
    {
    }

    /**
     * Compute the external force due to a Neumann type boundary condition.
     * This computes the following integral on a boundary element

     */
    virtual std::tuple<List const&, Vector> external_force(int const element,
                                                           double const load_factor) const override
    {
        auto const X = sf->project_to_plane(
            material_coordinates->initial_configuration(nodal_connectivity[element]));

        auto const h = interpolate_prescribed_load(load_factor);

        // Perform the computation of the external load vector
        auto const f_ext = sf->quadrature().integrate(Vector::Zero(X.cols()).eval(),
                                                      [&](auto const& femval, auto const& l) {
                                                          auto const & [ N, dN ] = femval;

                                                          auto const j = (X * dN).determinant();

                                                          return h * N * j;
                                                      });
        return {dof_list[element], f_ext};
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
                        std::shared_ptr<MaterialCoordinates>& material_coordinates,
                        Json::Value const& time_history,
                        Json::Value const& load_history,
                        int const dof_offset,
                        int const nodal_dofs)
        : Neumann(nodal_connectivity,
                  filter_dof_list(nodal_dofs, dof_offset, nodal_connectivity),
                  material_coordinates,
                  time_history,
                  load_history),
          sf(std::move(sf))
    {
    }

    std::tuple<List const&, Vector> external_force(int const element, double const load_factor) const
    {
        auto const X = material_coordinates->initial_configuration(nodal_connectivity[element]);

        auto const f = interpolate_prescribed_load(load_factor);

        // Perform the computation of the external load vector
        auto const f_ext = sf->quadrature().integrate(Vector::Zero(X.cols()).eval(),
                                                      [&](auto const& femval, auto const& l) {
                                                          auto const & [ N, dN ] = femval;

                                                          auto const j = (X * dN).determinant();

                                                          return f * N * j;
                                                      });
        return {dof_list.at(element), interpolate_prescribed_load(load_factor) * f_ext};
    }

protected:
    std::unique_ptr<VolumeInterpolation_Tp> sf;
};
}
