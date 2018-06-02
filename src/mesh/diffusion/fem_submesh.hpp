
#pragma once

#include "mesh/basic_submesh.hpp"

#include "constitutive/constitutive_model.hpp"
#include "constitutive/internal_variables.hpp"
#include "constitutive/internal_variables_alias.hpp"
#include "interpolations/shape_function.hpp"

#include <memory>

namespace neon
{
class material_coordinates;

namespace diffusion
{
/**
 * fem_submesh provides the element local routines for computing the system
 * components for a three-dimensional heat equation discretisation.
 */
class fem_submesh : public basic_submesh
{
public:
    using ValueCount = std::tuple<vector, vector>;

public:
    explicit fem_submesh(json const& material_data,
                         json const& mesh_data,
                         std::shared_ptr<material_coordinates>& mesh_coordinates,
                         basic_submesh const& submesh);

    /** @return list of global degrees of freedom for an element */
    [[nodiscard]] index_view local_dof_view(std::int64_t const element) const {
        return local_node_view(element);
    }

        /** @return The internal variable store */
        [[nodiscard]] auto const& internal_variables() const
    {
        return *variables;
    }

    void save_internal_variables(bool const have_converged);

    [[nodiscard]] auto dofs_per_node() const { return 1; }

    [[nodiscard]] auto const& shape_function() const { return *sf; }

    [[nodiscard]] auto const& constitutive() const { return *cm; }

    /**
     * Compute the stiffness (conductivity) matrix according to
     * \f{align*}{
     *     k_{ab} &= \int_{\Omega_e} \nabla N_a \kappa \nabla N_b d\Omega
     * \f}
     * where \f$ \kappa \f$ is the conductivity
     * @return DoFs and stiffness matrix
     */
    [[nodiscard]] std::pair<index_view, matrix> tangent_stiffness(std::int64_t const element) const;

    /**
     * Compute the consistent (full) mass matrix according to
     * \f{align*}{
     *     m_{ab} &= \int_{\Omega_e} N_a \rho c_p N_b d\Omega
     * \f}
     * where \f$ \rho \f$ is the density and \f$ c_p \f$ is the specific heat
     * @return DoFs and consistent mass matrix \sa diagonal_mass
     */
    [[nodiscard]] std::pair<index_view, matrix> consistent_mass(std::int64_t const element) const;

    /** @return Diagonal mass matrix using row sum technique \sa consistent_mass */
    [[nodiscard]] std::pair<index_view, vector> diagonal_mass(std::int64_t const element) const;

    /** Update the internal variables for the mesh group */
    void update_internal_variables(double const time_step_size);

    /**
     * Compute the local Jacobian matrix \f$ \bf{x}_\xi \f$
     * @param rhea Shape function gradients at quadrature point
     * @param configuration Configuration of the element (coordinates)
     */
    [[nodiscard]] matrix3 local_jacobian(matrix const& rhea, matrix const& configuration) const {
        return configuration * rhea;
    }

        [[nodiscard]] ValueCount
        nodal_averaged_variable(internal_variables_t::second const tensor_name) const;

    [[nodiscard]] ValueCount nodal_averaged_variable(internal_variables_t::scalar const scalar_name) const;

private:
    std::shared_ptr<material_coordinates> mesh_coordinates; /// Nodal coordinates

    std::unique_ptr<volume_interpolation> sf; /// Shape function

    variable_view view;
    std::shared_ptr<internal_variables_t> variables;

    std::unique_ptr<constitutive_model> cm; /// Constitutive model
};
}
}
