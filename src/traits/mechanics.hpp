
#pragma once

#include "numeric/dense_matrix.hpp"
#include "math/integral_form.hpp"

#include <array>
#include <type_traits>

namespace neon::mechanics
{
/// Specialisation of momentum equation for geometrical assumptions
enum class type {
    /// C1 beam theory (Hermite basis)
    euler_beam,
    /// C0 beam theory
    beam,
    shell,
    axisymmetric,
    plane_stress,
    plane_strain,
    continuum
};

enum class material { linear, nonlinear };

/// Deformation assumption
enum class deformation {
    /// Small displacements and rotations
    small,
    /// Finite displacements and rotations
    finite
};

enum class method {
    /// Linear perturbation is required and residual needs to be check manually
    linear,
    /// Nonlinear (incremental) solution of system
    incremental,
    /// Latin method
    latin
};

/// Frame
enum class lagrangian {
    /// Formulate integrals in reference configuration
    total,
    /// Formulate integrals in current configuration
    updated
};

/// Theoretical specialisations / approximations of Newton's laws.
enum class theory { plane_stress, plane_strain, axisymmetric, solid, shell, beam };

/// List of supported discretisations
enum class discretisation { linear, material_nonlinear, finite_strain, latin };

template <theory T, discretisation D, bool is_symmetric_ = true>
struct traits;

/// Trait specialisation for three dimensional beam theory
template <discretisation D, bool is_symmetric_>
struct traits<theory::beam, D, is_symmetric_>
{
    /// Theory
    static auto constexpr theory_type{theory::beam};
    /// Discretisation
    static auto constexpr discretisation_type{D};
    /// Spatial coordinates (x, y, z)
    static auto constexpr size{3};
    /// Is symmetric matrix
    static auto constexpr is_symmetric{is_symmetric_};
    /// Number of degrees of freedom per node
    static auto constexpr dofs_per_node{6};

    /// Three dimensional beams have three displacements followed by three rotations
    static std::array<int, dofs_per_node> constexpr dof_order{0, 1, 2, 3, 4, 5};

    using rank_two_tensor = matrix2;
    using rank_four_tensor = matrix2;
};

template <discretisation D, bool is_symmetric_>
struct traits<theory::plane_strain, D, is_symmetric_>
{
    static auto constexpr value_type{theory::plane_strain};
    // static auto constexpr type

    static auto constexpr size{2};
    static auto constexpr is_symmetric{is_symmetric_};
    static auto constexpr dofs_per_node{2};

    /// Plane stress mechanics has two unknowns (u_x, u_y)
    static std::array<int, dofs_per_node> constexpr dof_order{0, 1};

    using rank_two_tensor = matrix2;
    using rank_four_tensor = matrix3;
};

template <discretisation D, bool is_symmetric_>
struct traits<theory::plane_stress, D, is_symmetric_>
{
    static auto constexpr value{theory::plane_stress};
    static auto constexpr size{2};
    static auto constexpr is_symmetric{is_symmetric_};

    static auto constexpr dofs_per_node{2};

    /// Plane stress mechanics has two unknowns (u_x, u_y)
    static std::array<int, dofs_per_node> constexpr dof_order{0, 1};

    using rank_two_tensor = matrix2;
    using rank_four_tensor = matrix3;
};

template <discretisation D, bool is_symmetric_>
struct traits<theory::solid, D, is_symmetric_>
{
    static auto constexpr value{theory::solid};
    static auto constexpr size{3};
    static auto constexpr is_symmetric{is_symmetric_};

    static auto constexpr dofs_per_node{3};

    /// Solid mechanics has three unknowns (u_x, u_y, u_z)
    static std::array<int, dofs_per_node> constexpr dof_order{0, 1, 2};

    /// matrix type for the deformation and stress vector
    using rank_two_tensor = matrix3;
    /// matrix type for a fourth order tensor
    using rank_four_tensor = typename std::conditional<is_symmetric, matrix6, matrix9>::type;

    using bilinear_gradient = fem::integral<volume_interpolation, volume_quadrature, 1>;
    using bilinear = fem::integral<volume_interpolation, volume_quadrature, 0>;
};
}
