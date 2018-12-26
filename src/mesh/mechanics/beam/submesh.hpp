
#pragma once

#include "mesh/basic_submesh.hpp"

#include "mesh/material_coordinates.hpp"
#include "interpolations/shape_function.hpp"
#include "constitutive/internal_variables.hpp"

#include "constitutive/mechanics/beam/isotropic_linear.hpp"
#include "material/isotropic_elastic_property.hpp"
#include "math/view.hpp"

#include "geometry/profile.hpp"
#include "traits/mechanics.hpp"

#include <memory>
#include <utility>
#include <tuple>

namespace neon
{
class material_coordinates;
}

namespace neon::mechanics::beam
{
/// submesh is responsible for computing the linear stiffness matrices
/// for the three-dimensional small strain beam theory \cite Hughes2012.
class submesh : public basic_submesh //: public detail::submesh<beam::submesh>
{
public:
    // Type aliases
    using traits = mechanics::traits<theory::beam, discretisation::linear>;

    using internal_variable_type = internal_variables_t;

public:
    /// Constructor providing the material coordinates reference
    explicit submesh(json const& material_data,
                     json const& simulation_data,
                     json const& section_data,
                     std::shared_ptr<material_coordinates>& coordinates,
                     basic_submesh const& submesh);

    /// \return degrees of freedom and the linear element stiffness matrix
    [[nodiscard]] std::pair<index_view, matrix const&> tangent_stiffness(std::int32_t const element) const;

    /// Update the internal variables for the mesh group
    /// \sa update_deformation_measures()
    /// \sa update_Jacobian_determinants()
    /// \sa check_element_distortion()
    void update_internal_variables(double const time_step_size = 1.0);

    /// \return A view of degrees of freedom for an element
    [[nodiscard]] auto local_dof_view(std::int64_t const element) const
    {
        return dof_indices(Eigen::all, element);
    }

    /// \return The internal variable store
    [[nodiscard]] auto const& internal_variables() const { return *variables; }

    template <typename name_type>
    std::pair<vector, vector> nodal_averaged_variable(name_type) const
    {
        return std::make_pair(vector(), vector());
    }

protected:
    /**
     * Compute the bending contribution to the beam stiffness for an \p element
     * according to:
     * \f{align*}{
     *     \mathbf{k}^e_b &= \int_{0}^{h_e} \mathbf{B}^T_{b} \mathbf{D}_b \mathbf{B}_b \, dx_{3}^{e}
     * \f}
     * where \f$ \mathbf{B}_b \f$:
     * \f{align*}{
     * \mathbf{B}_b &= \begin{bmatrix} 0 & 0 & 0 & N_c^{'} & 0 & 0 \\
     *                                 0 & 0 & 0 & 0 & N_c^{'} & 0 \\
     * \end{bmatrix} \f}
     * and \f$ \mathbf{D}_b \f$ is:
     * \f{align*}{
     * \mathbf{D}_b &= \begin{bmatrix} E I_1 & 0 \\ 0 & E I_2 \end{bmatrix} \f}
     */
    [[nodiscard]] matrix const& bending_stiffness(matrix3x const& configuration,
                                                  std::int32_t const element) const;

    /**
     * Compute the shear contribution to the beam stiffness for an \p element
     * according to:
     * \f{align*}{
     *     \mathbf{k}^e_s &= \int_{0}^{h_e} \mathbf{B}^T_{s} \mathbf{D}_s \mathbf{B}_s \, dx_{3}^{e}
     * \f}
     * where \f$ \mathbf{B}_s \f$:
     * \f{align*}{
     * \mathbf{B}_s &= \begin{bmatrix} N_c^{'} & 0 & 0 & 0 & -N_c & 0 \\
     *                                 0 & N_c^{'} & 0 & N_c & 0 & 0\\
     * \end{bmatrix} \f}
     * and \f$ \mathbf{D}_s \f$ is:
     * \f{align*}{
     *     \mathbf{D}_s &= \begin{bmatrix} \mu A_1 & 0 \\ 0 & \mu A_1 \end{bmatrix}
     * \f}
     */
    [[nodiscard]] matrix const& shear_stiffness(matrix3x const& configuration,
                                                std::int32_t const element) const;

    /**
     * Compute the shear contribution to the beam stiffness for an \p element
     * according to:
     * \f{align*}{
     *     \mathbf{k}^e_s &= \int_{0}^{h_e} \mathbf{B}^T_{a} E A \mathbf{B}_{a} \, dx_{3}^{e}
     * \f}
     * where \f$ \mathbf{B}_s \f$:
     * \f{align*}{
     *     \mathbf{B}_a &= \begin{bmatrix} 0 & 0 & N_c^{'} & 0 & 0 & 0 \end{bmatrix}
     * \f}
     */
    [[nodiscard]] matrix const& axial_stiffness(matrix3x const& configuration,
                                                std::int32_t const element) const;

    /**
     * Compute the shear contribution to the beam stiffness for an \p element
     * according to:
     * \f{align*}{
     *     \mathbf{k}^e_t &= \int_{0}^{h_e} \mathbf{B}^T_{t} \mu J \mathbf{B}_{t} \, dx_{3}^{e}
     * \f}
     *
     * \f{align*}{\mathbf{B}_t &= \begin{bmatrix} 0 & 0 & 0 & 0 & 0 & N_c^{'} \end{bmatrix} \f}
     */
    [[nodiscard]] matrix const& torsional_stiffness(matrix3x const& configuration,
                                                    std::int32_t const element) const;

    /// Perform input checks and allocate normal and tangent vectors
    /// Throws \p std::domain_error if inputs do not exist or match expectations
    void allocate_normal_and_tangent(json const& profile_data);

protected:
    /// Line shape function
    std::unique_ptr<line_interpolation> sf;

    /// Coordinates of each mesh
    std::shared_ptr<material_coordinates> coordinates;

    /// View wrapper into internal variables
    stride_view<> view;
    /// Internal variables (geometry and constitutive)
    std::shared_ptr<internal_variable_type> variables;

    /// Constitutive model
    std::unique_ptr<isotropic_linear> cm;

    /// Local-global element indices map
    indices dof_indices;

    /// Element geometry profiles
    std::unique_ptr<geometry::profile> profile;
    /// Segment tangent vector
    vector3 tangent;
    /// Segment first normal
    vector3 normal;

    /// Rotation matrix for elements
    matrix12 rotation_matrix;
};
}
