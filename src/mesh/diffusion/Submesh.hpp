
#pragma once

#include "mesh/SubMesh.hpp"

#include "constitutive/ConstitutiveModel.hpp"
#include "constitutive/InternalVariables.hpp"
#include "interpolations/ShapeFunction.hpp"

#include <memory>

namespace neon
{
// Forward declarations
class ConstitutiveModel;
class MaterialCoordinates;

namespace diffusion
{
class femSubmesh : public SubMesh
{
public:
    using ValueCount = std::tuple<Vector, Vector>;

public:
    explicit femSubmesh(Json::Value const& material_data,
                        Json::Value const& simulation_data,
                        std::shared_ptr<MaterialCoordinates>& material_coordinates,
                        SubMesh const& submesh);

    /** @return list of global degrees of freedom for an element */
    List const& local_dof_list(int const element) const { return local_node_list(element); }

    /** @return The internal variable store */
    InternalVariables const& internal_variables() const { return variables; }

    void save_internal_variables(bool const have_converged);

    auto dofs_per_node() const { return 1; }

    auto const& shape_function() const { return *sf.get(); }

    auto const& constitutive() const { return *cm.get(); }

    /**
     * Compute the stiffness (conductivity) matrix according to
     * \f{align*}{
     *     k_{ab} &= \int_{\Omega_e} \nabla N_a \kappa \nabla N_b d\Omega
     * \f}
     * where \f$ \kappa \f$ is the conductivity
     * @return DoFs and stiffness matrix
     */
    std::tuple<List const&, Matrix> tangent_stiffness(int const element) const;

    /**
     * Compute the consistent (full) mass matrix according to
     * \f{align*}{
     *     m_{ab} &= \int_{\Omega_e} N_a \rho c_p N_b d\Omega
     * \f}
     * where \f$ \rho \f$ is the density and \f$ c_p \f$ is the specific heat
     * @return DoFs and consistent mass matrix \sa diagonal_mass
     */
    std::tuple<List const&, Matrix> consistent_mass(int const element) const;

    /** @return Diagonal mass matrix using row sum technique \sa consistent_mass */
    std::tuple<List const&, Vector> diagonal_mass(int const element) const;

    /** Update the internal variables for the mesh group */
    void update_internal_variables(double const time_step_size);

    /**
     * Compute the local Jacobian matrix \f$ \bf{x}_\xi \f$
     * @param rhea Shape function gradients at quadrature point
     * @param configuration Configuration of the element (coordinates)
     */
    Matrix3 local_jacobian(Matrix const& rhea, Matrix const& configuration) const
    {
        return configuration * rhea;
    }

    ValueCount nodal_averaged_variable(InternalVariables::Tensor const tensor_name) const;

    ValueCount nodal_averaged_variable(InternalVariables::Scalar const scalar_name) const;

protected:
    /** @return the index into the internal variable store */
    int offset(int const element, int const quadraturePoint) const;

private:
    std::shared_ptr<MaterialCoordinates> material_coordinates;

    std::unique_ptr<VolumeInterpolation> sf; //!< Shape function

    InternalVariables variables;

    std::unique_ptr<ConstitutiveModel> cm; //!< Constitutive model
};
}
}
