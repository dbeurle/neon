
#pragma once

#include "mesh/SubMesh.hpp"

#include "constitutive/ConstitutiveModel.hpp"
#include "constitutive/InternalVariables.hpp"
#include "interpolations/ShapeFunction.hpp"

#include <json/forwards.h>
#include <memory>
#include <unordered_map>

namespace neon
{
// Forward declarations
class ConstitutiveModel;

namespace solid
{
class MaterialCoordinates;

class femSubmesh : public SubMesh
{
public:
    using ValueCount = std::tuple<Vector, Vector>;

public:
    /** Constructor providing the material coordinates reference */
    explicit femSubmesh(Json::Value const& material_data,
                        Json::Value const& simulation_data,
                        std::shared_ptr<MaterialCoordinates>& material_coordinates,
                        SubMesh const& submesh);

    /** @return list of global degrees of freedom for an element */
    List const& local_dof_list(int element) const { return dof_list.at(element); }

    /** @return The internal variable store */
    InternalVariables const& internal_variables() const { return variables; }

    void save_internal_variables(bool const have_converged);

    auto dofs_per_node() const { return 3; }

    auto const& shape_function() const { return *sf.get(); }

    auto const& constitutive() const { return *cm.get(); }

    /** @return the tangent consistent stiffness matrix */
    std::tuple<List const&, Matrix> tangent_stiffness(int const element) const;

    /** @return the internal element force */
    std::tuple<List const&, Vector> internal_force(int const element) const;

    /** @return the consistent mass matrix \sa diagonal_mass */
    std::tuple<List const&, Matrix> consistent_mass(int const element) const;

    /** @return the consistent mass matrix \sa diagonal_mass */
    std::tuple<List const&, Vector> diagonal_mass(int const element) const;

    /** Update the internal variables for the mesh group
     *  Calls:
     *  \sa update_deformation_measures()
     *  \sa update_Jacobian_determinants()
     *  \sa check_element_distortion()
     */
    void update_internal_variables(double const time_step_size = 1.0);

    /**
     * Compute the local deformation gradient
     * \f{align*}{ F_{\xi} &= \bf{x}_\xi \f}
     * @param rhea Shape function gradients at quadrature point
     * @param configuration Configuration of the element (coordinates)
     */
    static Matrix3 local_deformation_gradient(Matrix const& rhea,
                                              Matrix const& configuration)
    {
        return configuration * rhea;
    }

    /**
     * Compute the deformation gradient, F, from the global to local mapping
     * \f{align*}{
     * F &= F_{\xi} \times (F^0_{\xi})^{-1}
     * \f}
     * @param rhea Shape function derivatives at the integration points
     * @param X Reference configuration (spatial coordinates, local nodes)
     * @param x Current configuration (spatial coordinates, local nodes)
     */
    Matrix3 deformation_gradient(Matrix const& rhea, Matrix const& X, Matrix const& x)
    {
        // Deformation gradient in the reference and current configuration
        return local_deformation_gradient(rhea, x)
               * local_deformation_gradient(rhea, X).inverse();
    }

    ValueCount nodal_averaged_variable(InternalVariables::Tensor const tensor_name) const;

    ValueCount nodal_averaged_variable(InternalVariables::Scalar const scalar_name) const;

protected:
    /** Update the strain measures defined by the constitutive model */
    void update_deformation_measures();

    /** Computes the Jacobian determinants.  Called by \sa update_deformation_measures()
     */
    void update_Jacobian_determinants();

    /** Check element Jacobians are acceptable.  Called by \sa
     * update_deformation_measures() */
    void check_element_distortion() const;

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
    Matrix geometric_tangent_stiffness(Matrix const& configuration, int element) const;

    /**
     * Compute the material tangent stiffness using the formula
     * \f{align*}{
     * k_{mat} &= I_{2x2} \int_{V} B_I^{T} \sigma B_{J} dV
     * \f}
     */
    Matrix material_tangent_stiffness(Matrix const& configuration, int element) const;

    /**
     * Compute the internal force vector using the formula
     * \f{align*}{
     * f_{int} &= \int_{V} B^{T} \sigma dV
     * \f}
     * @return the internal nodal force vector
     */
    Vector internal_nodal_force(Matrix const& configuration, int element) const;

    void allocate_dof_list(int const nodal_dofs);

    /** @return the index into the internal variable store */
    int offset(int element, int quadraturePoint) const;

private:
    std::shared_ptr<MaterialCoordinates> material_coordinates;

    std::unique_ptr<VolumeInterpolation> sf; //!< Shape function

    InternalVariables variables;

    std::unique_ptr<ConstitutiveModel> cm; //!< Constitutive model

    std::vector<List> dof_list; //!< Map for the local to global dofs
};
}
}
