
#pragma once

#include "mesh/Submesh.hpp"

#include "constitutive/ConstitutiveModel.hpp"
#include "constitutive/InternalVariables.hpp"
#include "interpolations/ShapeFunction.hpp"

#include <memory>

namespace neon
{
class MaterialCoordinates;

namespace diffusion
{
/**
 * femSubmesh provides the element local routines for computing the system
 * components for a three-dimensional heat equation discretisation.
 */
class femSubmesh : public Submesh
{
public:
    using ValueCount = std::tuple<vector, vector>;

public:
    explicit femSubmesh(Json::Value const& material_data,
                        Json::Value const& mesh_data,
                        std::shared_ptr<MaterialCoordinates>& material_coordinates,
                        Submesh const& submesh);

    /** @return list of global degrees of freedom for an element */
    [[nodiscard]] std::vector<int64> const& local_dof_list(int64 const element) const {
        return local_node_list(element);
    }

        /** @return The internal variable store */
        [[nodiscard]] InternalVariables const& internal_variables() const
    {
        return *variables;
    }

    void save_internal_variables(bool const have_converged);

    [[nodiscard]] auto dofs_per_node() const { return 1; }

    [[nodiscard]] auto const& shape_function() const { return *sf.get(); }

    [[nodiscard]] auto const& constitutive() const { return *cm.get(); }

    /**
     * Compute the stiffness (conductivity) matrix according to
     * \f{align*}{
     *     k_{ab} &= \int_{\Omega_e} \nabla N_a \kappa \nabla N_b d\Omega
     * \f}
     * where \f$ \kappa \f$ is the conductivity
     * @return DoFs and stiffness matrix
     */
    [[nodiscard]] std::pair<std::vector<int64> const&, matrix> tangent_stiffness(
        int64 const element) const;

    /**
     * Compute the consistent (full) mass matrix according to
     * \f{align*}{
     *     m_{ab} &= \int_{\Omega_e} N_a \rho c_p N_b d\Omega
     * \f}
     * where \f$ \rho \f$ is the density and \f$ c_p \f$ is the specific heat
     * @return DoFs and consistent mass matrix \sa diagonal_mass
     */
    [[nodiscard]] std::pair<std::vector<int64> const&, matrix> consistent_mass(int64 const element) const;

    /** @return Diagonal mass matrix using row sum technique \sa consistent_mass */
    [[nodiscard]] std::pair<std::vector<int64> const&, vector> diagonal_mass(int64 const element) const;

    /** Update the internal variables for the mesh group */
    void update_internal_variables(double const time_step_size);

    [[nodiscard]] ValueCount nodal_averaged_variable(InternalVariables::Tensor const tensor_name) const;

    [[nodiscard]] ValueCount nodal_averaged_variable(InternalVariables::Scalar const scalar_name) const;

protected:
    /** @return the index into the internal variable store */
    [[nodiscard]] int offset(int64 const element, int const quadraturePoint) const;

private:
    std::shared_ptr<MaterialCoordinates> material_coordinates;

    std::unique_ptr<VolumeInterpolation> sf; //!< Shape function

    std::shared_ptr<InternalVariables> variables;

    std::unique_ptr<ConstitutiveModel> cm; //!< Constitutive model
};
}
}
