
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
/// submesh provides the element local routines for computing the system
/// components for a three-dimensional heat equation discretisation.
class submesh : public basic_submesh
{
public:
    explicit submesh(json const& material_data,
                     json const& mesh_data,
                     std::shared_ptr<material_coordinates>& coordinates,
                     basic_submesh const& submesh);

    /// \return list of global degrees of freedom for an element
    [[nodiscard]] auto local_dof_view(std::int64_t const element) const -> index_view
    {
        return local_node_view(element);
    }

    /// \return The internal variable store
    [[nodiscard]] auto const& internal_variables() const { return *variables; }

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
    [[nodiscard]] auto tangent_stiffness(std::int64_t const element) const
        -> std::pair<index_view, matrix>;

    /**
     * Compute the consistent (full) mass matrix according to
     * \f{align*}{
     *     m_{ab} &= \int_{\Omega_e} N_a \rho c_p N_b d\Omega
     * \f}
     * where \f$ \rho \f$ is the density and \f$ c_p \f$ is the specific heat
     * @return DoFs and consistent mass matrix \sa diagonal_mass
     */
    [[nodiscard]] auto consistent_mass(std::int64_t const element) const
        -> std::pair<index_view, matrix>;

    /// \return Diagonal mass matrix using row sum technique \sa consistent_mass
    [[nodiscard]] auto diagonal_mass(std::int64_t const element) const
        -> std::pair<index_view, vector>;

    /// Update the internal variables for the mesh group
    void update_internal_variables(double const time_step_size);

    /// Compute the local Jacobian matrix \f$ \bf{x}_\xi \f$
    /// \param rhea Shape function gradients at quadrature point
    /// \param configuration Configuration of the element (coordinates)
    auto local_jacobian(matrix const& rhea, matrix const& configuration) const -> matrix3
    {
        return configuration * rhea;
    }

    ///
    [[nodiscard]] auto nodal_averaged_variable(variable::scalar const scalar_name) const
        -> std::pair<vector, vector>;

    ///
    [[nodiscard]] auto nodal_averaged_variable(variable::second const tensor_name) const
        -> std::pair<vector, vector>;

private:
    /// Nodal coordinates
    std::shared_ptr<material_coordinates> coordinates;
    /// Shape function
    std::unique_ptr<volume_interpolation> sf;

    variable_view view;
    std::shared_ptr<internal_variables_t> variables;

    /// Constitutive model
    std::unique_ptr<constitutive_model> cm;
};
}
}
