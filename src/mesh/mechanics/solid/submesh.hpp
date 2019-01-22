
#pragma once

/// @file

#include "mesh/basic_submesh.hpp"

#include "constitutive/constitutive_model.hpp"
#include "constitutive/internal_variables.hpp"
#include "math/integral_form.hpp"
#include "math/view.hpp"
#include "traits/mechanics.hpp"

#include <memory>

namespace neon
{
class material_coordinates;
class recovery;

namespace mechanics::solid
{
/// submesh provides the element local routines for computing the system
/// components for a three-dimensional continuum mechanics discretisation.
class submesh : public basic_submesh
{
public:
    using internal_variable_type = internal_variables_t;

    using traits = mechanics::traits<theory::solid, discretisation::finite_strain>;

    template <int N>
    using bilinear_type = fem::integral<volume_interpolation, volume_quadrature, N>;

public:
    explicit submesh(json const& material_data,
                     json const& mesh_data,
                     std::shared_ptr<material_coordinates>& coordinates,
                     basic_submesh const& submesh);

    ~submesh();

    submesh(submesh const&) = delete;

    submesh(submesh&&);

    submesh& operator=(submesh const&) = delete;

    submesh& operator=(submesh&&);

    /// \return view of degrees of freedom for an element
    [[nodiscard]] auto local_dof_view(std::int32_t const element) const
    {
        return dof_indices(Eigen::all, element);
    }

    /// \return internal variables at the quadrature point
    [[nodiscard]] auto const& internal_variables() const { return *variables; }

    void save_internal_variables(bool const have_converged);

    /// \return number of degrees of freedom per node
    [[nodiscard]] auto dofs_per_node() const noexcept { return traits::dofs_per_node; }

    /// \return underlying constitutive model
    [[nodiscard]] auto const& constitutive() const { return *cm; }

    /// \return tangent consistent stiffness matrix
    [[nodiscard]] matrix const& tangent_stiffness(std::int32_t const element) const;

    /**
     * Compute the internal force vector using the formula
     * \f{align*}{
     * f_{int} &= \int_{V} B^{T} \sigma dV
     * \f}
     * \return internal element force
     */
    [[nodiscard]] vector const& internal_force(std::int32_t const element) const;

    /// \return consistent mass matrix \sa diagonal_mass
    [[nodiscard]] matrix const& consistent_mass(std::int32_t const element) const;

    /// \return consistent mass matrix \sa diagonal_mass
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

    /// Compute the Jacobian determinants and check if negative
    void update_jacobian_determinants();

    /**
     * Compute the geometric stiffness matrix for the solid element. The
     * expression to be evaluated through numerical integration is:
       \f{align*}{
        k_{geo} &= \sum_l^{L} B(\xi_l, \eta_l, \zeta_l)^T \sigma(l) B(\xi_l,
     \eta_l, \zeta_l) w(l)
       \f}
     * Where B is the gradient operator in the finite element discretization
     */
    [[nodiscard]] matrix const& geometric_tangent_stiffness(matrix3x const& configuration,
                                                            std::int32_t const element) const;

    /**
     * Compute the material tangent stiffness using the formula
     * \f{align*}{
     * k_{mat} &= I_{2x2} \int_{V} B_I^{T} \sigma B_{J} dV
     * \f}
     */
    [[nodiscard]] matrix const& material_tangent_stiffness(matrix3x const& configuration,
                                                           std::int32_t const element) const;

protected:
    std::shared_ptr<material_coordinates> coordinates;

    /// Gradient bilinear form for stiffness matrix
    bilinear_type<0> bilinear_gradient;
    /// Bilinear form for mass matrix
    bilinear_type<0> bilinear;

    stride_view<> view;
    std::shared_ptr<internal_variable_type> variables;
    /// Constitutive model
    std::unique_ptr<constitutive_model> cm;
    /// Map for the local to global dofs
    indices dof_indices;

    std::unique_ptr<recovery> patch_recovery;
};
}
}
