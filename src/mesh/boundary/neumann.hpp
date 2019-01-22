
#pragma once

#include "boundary_condition.hpp"

#include "io/json_forward.hpp"
#include "math/jacobian_determinant.hpp"
#include "mesh/element_topology.hpp"
#include "mesh/material_coordinates.hpp"
#include "numeric/index_types.hpp"

#include <memory>

namespace neon
{
/// neumann is a base class for Neumann (derivative) type boundary conditions.
/// This includes the nodal connectivities and degrees of freedom lists.  Derived
/// classes must implement shape functions and the appropriate finite element
/// approximation for the given problem
class neumann : public boundary_condition
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

    [[nodiscard]] auto local_dof_view(std::int64_t const element) const noexcept
    {
        return dof_indices(Eigen::all, element);
    }

    [[nodiscard]] auto local_node_view(std::int64_t const element) const noexcept
    {
        return node_indices(Eigen::all, element);
    }

protected:
    /// Indices for the nodal coordinates
    indices node_indices;
    /// Indices for the degrees of freedom
    indices dof_indices;
    /// Coordinates for the boundary element group
    std::shared_ptr<material_coordinates> coordinates;
};

/// constant_neumann is a specialisation of a \p neumann boundary condition that
/// computes a surface integral for scalar loads
template <typename LinearFormType>
class constant_neumann : public neumann
{
public:
    /// Type alias for the integration form
    using linear_form_type = LinearFormType;

public:
    explicit constant_neumann(indices node_indices,
                              indices dof_indices,
                              std::shared_ptr<material_coordinates>& coordinates,
                              json const& time_history,
                              json const& load_history,
                              element_topology const topology,
                              json const& element_options)
        : neumann(node_indices, dof_indices, coordinates, time_history, load_history),
          linear_form(topology, element_options)
    {
    }

    explicit constant_neumann(indices node_indices,
                              indices dof_indices,
                              std::shared_ptr<material_coordinates>& coordinates,
                              json const& boundary,
                              std::string const& name,
                              double const generate_time_step,
                              element_topology const topology,
                              json const& element_options)
        : neumann(node_indices, dof_indices, coordinates, boundary, name, generate_time_step),
          linear_form(topology, element_options)
    {
    }

    virtual ~constant_neumann() = default;

    constant_neumann(constant_neumann&&) = default;

    constant_neumann& operator=(constant_neumann const&) = default;

    /// Compute the external force contribution due to a neumann type boundary
    /// condition. This computes the following integral on a boundary element
    /// \param element Element to compute external force
    /// \param load_factor Load factor to interpolate load
    /// \return Load vector for assembly
    virtual auto external_force(std::int64_t const element, double const load_factor) const
        -> vector const&
    {
        thread_local vector f_ext;

        auto const node_view = local_node_view(element);

        auto const X = coordinates->initial_configuration(node_view);

        f_ext.resize(node_view.size());

        // Perform the computation of the external load vector
        linear_form.integrate(f_ext.setZero(), [&](auto const& value, auto) -> vector {
            auto const& [N, dN] = value;

            return N * jacobian_determinant(X * dN);
        });

        f_ext *= interpolate_prescribed_load(load_factor);

        return f_ext;
    }

protected:
    /// Linear form for the integral evaluation
    linear_form_type linear_form;
};
}
