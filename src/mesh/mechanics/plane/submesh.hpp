
#pragma once

/// @file

#include "mesh/mechanics/submesh.hpp"

#include "constitutive/constitutive_model.hpp"
#include "constitutive/internal_variables.hpp"
#include "mesh/projection/recovery.hpp"
#include "math/view.hpp"
#include "math/integral_form.hpp"
#include "traits/mechanics.hpp"

#include <memory>
#include <utility>

namespace neon
{
class material_coordinates;

namespace mechanics::plane
{
/// submesh provides the element local routines for computing the system
/// components for a generalised plane strain/stress mechanics discretisation.
/// This class conforms to the CRTP interface \sa detail::submesh
class submesh : public detail::submesh<plane::submesh, plane::internal_variables_t>
{
public:
    using traits = mechanics::traits<theory::plane_strain, discretisation::linear>;

public:
    /// Constructor providing the material coordinates reference
    explicit submesh(json const& material_data,
                     json const& simulation_data,
                     std::shared_ptr<material_coordinates>& coordinates,
                     basic_submesh const& submesh);

    submesh(submesh&&) = default;

    [[nodiscard]] auto local_dof_view(std::int32_t const element) const noexcept
    {
        return dof_list(Eigen::all, element);
    }

    [[nodiscard]] auto const& internal_variables() const { return *variables; }

    void save_internal_variables(bool const have_converged);

    [[nodiscard]] auto dofs_per_node() const noexcept { return traits::dofs_per_node; }

    [[nodiscard]] auto const& constitutive() const { return *cm; }

    /// \return the tangent consistent stiffness matrix
    [[nodiscard]] matrix const& tangent_stiffness(std::int32_t const element) const;

    /// \return the internal element force
    [[nodiscard]] vector const& internal_force(std::int32_t const element) const;

    /// \return the consistent mass matrix \sa diagonal_mass
    [[nodiscard]] matrix const& consistent_mass(std::int32_t const element) const;

    /// \return the consistent mass matrix \sa diagonal_mass
    [[nodiscard]] vector const& diagonal_mass(std::int32_t const element) const;

    /// Update the internal variables for the mesh group
    /// \sa update_deformation_measures()
    /// \sa update_jacobian_determinants()
    /// \sa check_element_distortion()
    void update_internal_variables(double const time_step_size = 1.0);

    [[nodiscard]] std::pair<vector, vector> nodal_averaged_variable(
        variable::second const tensor_name) const;

    [[nodiscard]] std::pair<vector, vector> nodal_averaged_variable(
        variable::scalar const scalar_name) const;

protected:
    /// Update the strain measures defined by the constitutive model
    void update_deformation_measures();

    /// Computes the Jacobian determinants and check if negative
    void update_jacobian_determinants();

    /**
     * Compute the geometric stiffness matrix for in the computation solid
     * mechanics element. The expression to be evaluated through numerical
     * integration is the following:
       \f{align*}{
        k_{geo} &= \sum_l^{L} B(\xi_l, \eta_l, \zeta_l)^T \sigma(l) B(\xi_l,
     \eta_l, \zeta_l) w(l)
       \f}
     * Where B is the gradient operator in the finite element discretization
     */
    [[nodiscard]] matrix const& geometric_tangent_stiffness(matrix2x const& configuration,
                                                            std::int32_t const element) const;

    /**
     * Compute the material tangent stiffness using the formula
     * \f{align*}{
     * k_{mat} &= I_{2x2} \int_{V} B_I^{T} \sigma B_{J} dV
     * \f}
     */
    [[nodiscard]] matrix const& material_tangent_stiffness(matrix2x const& configuration,
                                                           std::int32_t const element) const;

    /**
     * Compute the internal force vector using the formula
     * \f{align*}{
     * f_{int} &= \int_{V} B^{T} \sigma dV
     * \f}
     * @return the internal nodal force vector
     */
    [[nodiscard]] vector const& internal_nodal_force(matrix2x const& configuration,
                                                     std::int32_t const element) const;

private:
    std::shared_ptr<material_coordinates> coordinates;

    /// Bilinear gradient form for stiffness matrix and internal force
    fem::integral<surface_interpolation, surface_quadrature, 1> bilinear_gradient;
    /// Bilinear form for mass matrix
    fem::integral<surface_interpolation, surface_quadrature, 0> bilinear;

    stride_view<> view;
    std::shared_ptr<internal_variables_t> variables;

    /// Constitutive model
    std::unique_ptr<constitutive_model> cm;

    /// Map for the local element to process indices
    indices dof_list;

    std::unique_ptr<local_extrapolation> patch_recovery = nullptr;
};
}
}
