
#pragma once

#include "vector_contribution.hpp"

#include "math/jacobian_determinant.hpp"
#include "mesh/material_coordinates.hpp"

#include <memory>

namespace neon
{
/// neumann is a base class for Neumann (derivative) type boundary conditions.
/// This includes the nodal connectivities and degrees of freedom lists.  Derived
/// classes must implement shape functions and the appropriate finite element
/// approximation for the given problem
class neumann : public vector_contribution
{
public:
    /// \param node_indices Nodal list
    /// \param dof_indices Degree of freedom mapping
    /// \param material_coordinates Nodal coordinates
    /// \param times Input times
    /// \param times Input boundary values for each time
    explicit neumann(indices node_indices,
                     indices dof_indices,
                     std::shared_ptr<material_coordinates>& material_coordinates,
                     json const& times,
                     json const& loads);

    /// \param node_indices Nodal list
    /// \param dof_indices Degree of freedom mapping
    /// \param material_coordinates Nodal coordinates
    /// \param boundary Boundary information
    /// \param
    explicit neumann(indices node_indices,
                     indices dof_indices,
                     std::shared_ptr<material_coordinates>& material_coordinates,
                     json const& boundary,
                     std::string const& name,
                     double const generate_time_step);

    [[nodiscard]] auto elements() const noexcept { return node_indices.cols(); }

protected:
    /// Indices for the nodal coordinates
    indices node_indices;
    /// Indices for the degrees of freedom
    indices dof_indices;
    /// Coordinates for the boundary element group
    std::shared_ptr<material_coordinates> coordinates;
};

/// surface_load is a specialisation of a \p neumann boundary condition that
/// computes a surface integral for scalar loads
template <typename surface_interpolation_type>
class surface_load : public neumann
{
public:
    /// Type alias for the surface interpolation type
    using surface_interpolation = surface_interpolation_type;

public:
    explicit surface_load(std::unique_ptr<surface_interpolation>&& sf,
                          indices node_indices,
                          indices dof_indices,
                          std::shared_ptr<material_coordinates>& coordinates,
                          json const& time_history,
                          json const& load_history)
        : neumann(node_indices, dof_indices, coordinates, time_history, load_history),
          sf(std::move(sf))
    {
    }

    explicit surface_load(std::unique_ptr<surface_interpolation>&& sf,
                          indices node_indices,
                          indices dof_indices,
                          std::shared_ptr<material_coordinates>& coordinates,
                          json const& boundary,
                          std::string const& name,
                          double const generate_time_step)
        : neumann(node_indices, dof_indices, coordinates, boundary, name, generate_time_step),
          sf(std::move(sf))
    {
    }

    virtual ~surface_load() = default;

    surface_load(surface_load&& other) = default;
    surface_load& operator=(surface_load const&) = default;

    /// Compute the external force due to a neumann type boundary condition.
    /// This computes the following integral on a boundary element
    /// \param element Surface element to compute external force on
    /// \param load_factor Load factor to interpolate load with
    /// \return Dof list and a vector for assembly
    virtual std::pair<index_view, vector> external_force(std::int64_t const element,
                                                         double const load_factor) const override
    {
        auto const node_view = node_indices(Eigen::placeholders::all, element);

        auto const X = coordinates->initial_configuration(node_view);

        // Perform the computation of the external load vector
        auto const f_ext = sf->quadrature().integrate(vector::Zero(X.cols()).eval(),
                                                      [&](auto const& femval, auto) -> vector {
                                                          auto const& [N, dN] = femval;

                                                          return N * jacobian_determinant(X * dN);
                                                      });

        return {dof_indices(Eigen::placeholders::all, element),
                interpolate_prescribed_load(load_factor) * f_ext};
    }

protected:
    /// Shape function for surface interpolation
    std::unique_ptr<surface_interpolation> sf;
};

/// volume_load provides an interface for computing the contribution to the
/// load vector from a volume load boundary.
template <typename volume_interpolation_type>
class volume_load : public neumann
{
public:
    /// Type alias for the volume interpolation type
    using volume_interpolation = volume_interpolation_type;

public:
    explicit volume_load(std::unique_ptr<volume_interpolation>&& sf,
                         indices node_indices,
                         indices dof_indices,
                         std::shared_ptr<material_coordinates>& coordinates,
                         json const& time_history,
                         json const& load_history)
        : neumann(node_indices, dof_indices, coordinates, time_history, load_history),
          sf(std::move(sf))
    {
    }

    explicit volume_load(std::unique_ptr<volume_interpolation>&& sf,
                         indices node_indices,
                         indices dof_indices,
                         std::shared_ptr<material_coordinates>& coordinates,
                         json const& boundary,
                         std::string const& name,
                         double const generate_time_step)
        : neumann(node_indices, dof_indices, coordinates, boundary, name, generate_time_step),
          sf(std::move(sf))
    {
    }

    virtual ~volume_load() = default;

    volume_load(volume_load&& other) = default;

    volume_load& operator=(volume_load const&) = default;

    virtual std::pair<index_view, vector> external_force(std::int64_t const element,
                                                         double const load_factor) const override
    {
        auto const node_view = node_indices(Eigen::placeholders::all, element);

        auto const X = coordinates->initial_configuration(node_view);

        // Perform the computation of the external load vector
        auto const f_ext = sf->quadrature().integrate(vector::Zero(X.cols()).eval(),
                                                      [&](auto const& femval, auto) -> vector {
                                                          auto const& [N, dN] = femval;

                                                          return N * jacobian_determinant(X * dN);
                                                      });

        return {dof_indices(Eigen::placeholders::all, element),
                interpolate_prescribed_load(load_factor) * f_ext};
    }

protected:
    std::unique_ptr<volume_interpolation> sf;
};
}
