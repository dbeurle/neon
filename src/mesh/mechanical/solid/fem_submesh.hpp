
#pragma once

#include "mesh/basic_submesh.hpp"

#include "constitutive/constitutive_model.hpp"
#include "constitutive/internal_variables.hpp"
#include "interpolations/shape_function.hpp"
#include "traits/mechanics.hpp"

#include <memory>

namespace neon
{
class material_coordinates;

namespace mechanical::solid
{
/// fem_submesh provides the element local routines for computing the system
/// components for a three-dimensional continuum mechanics discretisation.
class fem_submesh : public basic_submesh
{
public:
    using ValueCount = std::pair<vector, vector>;

    using internal_variable_type = internal_variables_t;

    using trait_type = traits<theory::solid, discretisation::finite_strain>;

public:
    /** Constructor providing the material coordinates reference */
    explicit fem_submesh(json const& material_data,
                         json const& mesh_data,
                         std::shared_ptr<material_coordinates>& material_coordinates,
                         basic_submesh const& submesh);

    /// \return list of global degrees of freedom for an element
    [[nodiscard]] auto const local_dof_view(std::int64_t const element) const
    {
        return dof_list(Eigen::placeholders::all, element);
    }

    /// \return The internal variable store
    [[nodiscard]] auto const& internal_variables() const { return *variables; }

    void save_internal_variables(bool const have_converged);

    [[nodiscard]] auto dofs_per_node() const noexcept { return 3; }

    [[nodiscard]] auto const& shape_function() const { return *sf; }

    [[nodiscard]] auto const& constitutive() const { return *cm; }

    /** @return the tangent consistent stiffness matrix */
    [[nodiscard]] std::pair<index_view, matrix> tangent_stiffness(std::int32_t const element) const;

    /** @return the internal element force */
    [[nodiscard]] std::pair<index_view, vector> internal_force(std::int32_t const element) const;

    /** @return the consistent mass matrix \sa diagonal_mass */
    [[nodiscard]] std::pair<index_view, matrix> consistent_mass(std::int32_t const element) const;

    /** @return the consistent mass matrix \sa diagonal_mass */
    [[nodiscard]] std::pair<index_view, vector> diagonal_mass(std::int32_t const element) const;

    /// Update the internal variables for the mesh group
    /// \sa update_deformation_measures()
    /// \sa update_Jacobian_determinants()
    /// \sa check_element_distortion()
    void update_internal_variables(double const time_step_size = 1.0);

    [[nodiscard]] ValueCount nodal_averaged_variable(internal_variables_t::Tensor const tensor_name) const;

    [[nodiscard]] ValueCount nodal_averaged_variable(internal_variables_t::scalar const scalar_name) const;

protected:
    /** Update the strain measures defined by the constitutive model */
    void update_deformation_measures();

    /** Computes the Jacobian determinants and check if negative
     */
    void update_Jacobian_determinants();

    /**
     * Compute the geometric stiffness matrix for the solid element. The
     * expression to be evaluated through numerical integration is:
       \f{align*}{
        k_{geo} &= \sum_l^{L} B(\xi_l, \eta_l, \zeta_l)^T \sigma(l) B(\xi_l,
     \eta_l, \zeta_l) w(l)
       \f}
     * Where B is the gradient operator in the finite element discretization
     */
    [[nodiscard]] matrix geometric_tangent_stiffness(matrix3x const& configuration,
                                                     std::int32_t const element) const;

    /**
     * Compute the material tangent stiffness using the formula
     * \f{align*}{
     * k_{mat} &= I_{2x2} \int_{V} B_I^{T} \sigma B_{J} dV
     * \f}
     */
    [[nodiscard]] matrix material_tangent_stiffness(matrix3x const& configuration,
                                                    std::int32_t const element) const;

    /**
     * Compute the internal force vector using the formula
     * \f{align*}{
     * f_{int} &= \int_{V} B^{T} \sigma dV
     * \f}
     * @return the internal nodal force vector
     */
    [[nodiscard]] vector internal_nodal_force(matrix3x const& configuration,
                                              std::int32_t const element) const;

private:
    std::shared_ptr<material_coordinates> mesh_coordinates;

    std::unique_ptr<volume_interpolation> sf; //!< Shape function

    variable_view view;
    std::shared_ptr<internal_variables_t> variables;

    std::unique_ptr<constitutive_model> cm; //!< Constitutive model

    indices dof_list; //!< Map for the local to global dofs
};
}
}
