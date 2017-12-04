
#pragma once

#include "mesh/Submesh.hpp"

#include "constitutive/ConstitutiveModel.hpp"
#include "constitutive/InternalVariables.hpp"
#include "interpolations/ShapeFunction.hpp"

#include <memory>
#include <unordered_map>

namespace neon
{
class MaterialCoordinates;

namespace mech::solid
{
/**
 * femSubmesh provides the element local routines for computing the system
 * components for a three-dimensional continuum mechanics discretisation.
 */
class femSubmesh : public Submesh
{
public:
    using ValueCount = std::tuple<Vector, Vector>;

public:
    /** Constructor providing the material coordinates reference */
    explicit femSubmesh(Json::Value const& material_data,
                        Json::Value const& simulation_data,
                        std::shared_ptr<MaterialCoordinates>& material_coordinates,
                        Submesh const& submesh);

    /** @return list of global degrees of freedom for an element */
    [[nodiscard]] List const& local_dof_list(int const element) const
    {
        return dof_list.at(element);
    }

    /** @return The internal variable store */
    [[nodiscard]] InternalVariables const& internal_variables() const { return variables; }

    void save_internal_variables(bool const have_converged);

    [[nodiscard]] auto dofs_per_node() const { return 3; }

    [[nodiscard]] auto const& shape_function() const { return *sf.get(); }

    [[nodiscard]] auto const& constitutive() const { return *cm.get(); }

    /** @return the tangent consistent stiffness matrix */
    [[nodiscard]] std::tuple<List const&, Matrix> tangent_stiffness(int const element) const;

    /** @return the internal element force */
    [[nodiscard]] std::tuple<List const&, Vector> internal_force(int const element) const;

    /** @return the consistent mass matrix \sa diagonal_mass */
    [[nodiscard]] std::tuple<List const&, Matrix> consistent_mass(int const element) const;

    /** @return the consistent mass matrix \sa diagonal_mass */
    [[nodiscard]] std::tuple<List const&, Vector> diagonal_mass(int const element) const;

    /** Update the internal variables for the mesh group
     *  \sa update_deformation_measures()
     *  \sa update_Jacobian_determinants()
     *  \sa check_element_distortion()
     */
    void update_internal_variables(double const time_step_size = 1.0);

    [[nodiscard]] ValueCount nodal_averaged_variable(InternalVariables::Tensor const tensor_name) const;

    [[nodiscard]] ValueCount nodal_averaged_variable(InternalVariables::Scalar const scalar_name) const;

protected:
    /** Update the strain measures defined by the constitutive model */
    void update_deformation_measures();

    /** Computes the Jacobian determinants and check if negative
     */
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
    [[nodiscard]] Matrix geometric_tangent_stiffness(Matrix const& configuration, int element) const;

    /**
     * Compute the material tangent stiffness using the formula
     * \f{align*}{
     * k_{mat} &= I_{2x2} \int_{V} B_I^{T} \sigma B_{J} dV
     * \f}
     */
    [[nodiscard]] Matrix material_tangent_stiffness(Matrix const& configuration, int element) const;

    /**
     * Compute the internal force vector using the formula
     * \f{align*}{
     * f_{int} &= \int_{V} B^{T} \sigma dV
     * \f}
     * @return the internal nodal force vector
     */
    [[nodiscard]] Vector internal_nodal_force(Matrix const& configuration, int element) const;

    /** @return the index into the internal variable store */
    [[nodiscard]] int offset(int const element, int const quadraturePoint) const;

private:
    std::shared_ptr<MaterialCoordinates> material_coordinates;

    std::unique_ptr<VolumeInterpolation> sf; //!< Shape function

    InternalVariables variables;

    std::unique_ptr<ConstitutiveModel> cm; //!< Constitutive model

    std::vector<List> dof_list; //!< Map for the local to global dofs
};
}
}
