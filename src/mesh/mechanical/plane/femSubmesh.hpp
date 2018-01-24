
#pragma once

#include "mesh/mechanical/femSubmesh.hpp"

#include "constitutive/ConstitutiveModel.hpp"
#include "constitutive/InternalVariables.hpp"
#include "interpolations/shape_function.hpp"

#include <memory>

namespace neon
{
class MaterialCoordinates;

namespace mechanical::plane
{
/**
 * femSubmesh provides the element local routines for computing the system
 * components for a generalised plane strain/stress mechanics discretisation.
 * This class conforms to the CRTP interface \sa detail::femSubmesh
 */
class femSubmesh : public detail::femSubmesh<plane::femSubmesh, plane::InternalVariables>
{
public:
    /** Constructor providing the material coordinates reference */
    explicit femSubmesh(json const& material_data,
                        json const& simulation_data,
                        std::shared_ptr<MaterialCoordinates>& material_coordinates,
                        Submesh const& submesh);

    femSubmesh(femSubmesh&&) = default;

    [[nodiscard]] local_indices const& local_dof_view(int32 const element) const {
        return dof_list.at(element);
    }

        [[nodiscard]] InternalVariables const& internal_variables() const
    {
        return *variables;
    }

    void save_internal_variables(bool const have_converged);

    [[nodiscard]] auto dofs_per_node() const { return 2; }

    [[nodiscard]] auto const& shape_function() const { return *sf.get(); }

    [[nodiscard]] auto const& constitutive() const { return *cm.get(); }

    /** @return the tangent consistent stiffness matrix */
    [[nodiscard]] std::pair<local_indices const&, matrix> tangent_stiffness(int32 const element) const;

    /** @return the internal element force */
    [[nodiscard]] std::pair<local_indices const&, vector> internal_force(int32 const element) const;

    /** @return the consistent mass matrix \sa diagonal_mass */
    [[nodiscard]] std::pair<local_indices const&, matrix> consistent_mass(int32 const element) const;

    /** @return the consistent mass matrix \sa diagonal_mass */
    [[nodiscard]] std::pair<local_indices const&, vector> diagonal_mass(int32 const element) const;

    /** Update the internal variables for the mesh group
     *  \sa update_deformation_measures()
     *  \sa update_Jacobian_determinants()
     *  \sa check_element_distortion()
     */
    void update_internal_variables(double const time_step_size = 1.0);

    [[nodiscard]] std::pair<vector, vector> nodal_averaged_variable(
        InternalVariables::Tensor const tensor_name) const;

    [[nodiscard]] std::pair<vector, vector> nodal_averaged_variable(
        InternalVariables::Scalar const scalar_name) const;

protected:
    /** Update the strain measures defined by the constitutive model */
    void update_deformation_measures();

    /** Computes the Jacobian determinants and check if negative */
    void update_Jacobian_determinants();

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
    [[nodiscard]] matrix geometric_tangent_stiffness(matrix2x const& configuration,
                                                     int32 const element) const;

    /**
     * Compute the material tangent stiffness using the formula
     * \f{align*}{
     * k_{mat} &= I_{2x2} \int_{V} B_I^{T} \sigma B_{J} dV
     * \f}
     */
    [[nodiscard]] matrix material_tangent_stiffness(matrix2x const& configuration,
                                                    int32 const element) const;

    /**
     * Compute the internal force vector using the formula
     * \f{align*}{
     * f_{int} &= \int_{V} B^{T} \sigma dV
     * \f}
     * @return the internal nodal force vector
     */
    [[nodiscard]] vector internal_nodal_force(matrix2x const& configuration,
                                              int32 const element) const;

private:
    std::shared_ptr<MaterialCoordinates> material_coordinates;

    std::unique_ptr<surface_interpolation> sf; //!< Shape function

    variable_view view;
    std::shared_ptr<InternalVariables> variables;

    std::unique_ptr<ConstitutiveModel> cm; //!< Constitutive model

    std::vector<local_indices> dof_list; //!< Map for the local element to process indices
};
}
}
