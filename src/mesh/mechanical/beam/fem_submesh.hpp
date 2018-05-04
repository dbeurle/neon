
#pragma once

#include "interpolations/shape_function.hpp"
#include "constitutive/internal_variables.hpp"
#include "constitutive/mechanical/beam/isotropic_linear.hpp"
#include "material/isotropic_elastic_property.hpp"
#include "mesh/basic_submesh.hpp"
#include "mesh/material_coordinates.hpp"
#include "mesh/nodal_variables.hpp"
#include "numeric/dense_matrix.hpp"
#include "numeric/index_types.hpp"
#include "geometry/profile.hpp"
#include "traits/mechanics.hpp"

#include <memory>
#include <utility>
#include <tuple>

namespace neon
{
class material_coordinates;
}

namespace neon::mechanical::beam
{
/// fem_submesh is responsible for computing the linear stiffness matrices
/// for the three-dimensional small strain beam theory \cite Hughes2012.
class fem_submesh : public basic_submesh //: public detail::fem_submesh<beam::fem_submesh>
{
public:
    // Type aliases
    using traits = mechanical::traits<theory::beam, discretisation::linear>;

    using internal_variable_type = internal_variables<traits::rank_two_tensor::RowsAtCompileTime,
                                                      traits::rank_four_tensor::RowsAtCompileTime>;

public:
    /// Constructor providing the material coordinates reference
    explicit fem_submesh(json const& material_data,
                         json const& simulation_data,
                         std::shared_ptr<material_coordinates>& coordinates,
                         basic_submesh const& submesh);

    /// \return degrees of freedom and the linear element stiffness matrix
    [[nodiscard]] std::pair<index_view, matrix> tangent_stiffness(std::int32_t const element) const;

    /// Update the internal variables for the mesh group
    /// \sa update_deformation_measures()
    /// \sa update_Jacobian_determinants()
    /// \sa check_element_distortion()
    void update_internal_variables(double const time_step_size = 1.0);

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
    matrix const& bending_stiffness(matrix const& configuration, std::int32_t const element) const;

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
    matrix const& shear_stiffness(matrix const& configuration, std::int32_t const element) const;

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
    matrix const& axial_stiffness(matrix const& configuration, std::int32_t const element) const;

    /**
     * Compute the shear contribution to the beam stiffness for an \p element
     * according to:
     * \f{align*}{
     *     \mathbf{k}^e_t &= \int_{0}^{h_e} \mathbf{B}^T_{t} \mu J \mathbf{B}_{t} \, dx_{3}^{e}
     * \f}
     *
     * \f{align*}{\mathbf{B}_t &= \begin{bmatrix} 0 & 0 & 0 & 0 & 0 & N_c^{'} \end{bmatrix} \f}
     */
    matrix const& torsional_stiffness(matrix const& configuration, std::int32_t const element) const;

protected:
    /// Material properties
    isotropic_elastic_property material;

    /// Line shape function (linear, quadratic, cubic)
    std::unique_ptr<line_interpolation> sf;

    /// Coordinates of each mesh
    std::shared_ptr<material_coordinates> coordinates;

    /// u1, u2, u3, theta1, theta2, theta3
    nodal_variables<traits::dofs_per_node> displacement_rotation;

    /// Element profiles
    std::vector<std::unique_ptr<geometry::profile>> profiles;

    /// View wrapper into internal variables
    variable_view view;
    /// Internal variables (geometry and constitutive)
    std::shared_ptr<internal_variable_type> variables;

    std::unique_ptr<isotropic_linear> cm;

    /// Local-global element indices map
    indices dof_list;

    /// Element tangent vectors
    std::vector<vector3> tangents;
    /// Element orientiations
    std::vector<vector3> orientations;
};
}
