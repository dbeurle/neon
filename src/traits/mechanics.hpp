
#pragma once

#include <cstdint>

namespace neon::mechanics
{
enum class type : uint8_t { plane_stress, plane_strain, axisymmetric, solid, plate, shell };

template <neon::mechanics::type T, bool is_symmetric_>
struct traits;

template <>
struct traits<type::plane_strain, bool is_symmetric_ = true>
{
    static auto constexpr value = type::plane_strain;

    static auto constexpr size = 2;

    static auto constexpr is_symmetric = is_symmetric_;

    using constitutive_type = constitutive_model<value>;

    using mesh_type = fem_mesh<value>;

    using rank_two_tensor = matrix2;
    using rank_four_tensor = matrix3;
};

template <>
struct traits<type::plane_stress, bool is_symmetric_ = true>
{
    static auto constexpr value = type::plane_stress;

    static auto constexpr size = 2;

    static auto constexpr is_symmetric = is_symmetric_;

    using constitutive_type = constitutive_model<value>;

    using mesh_type = fem_mesh<value>;

    using rank_two_tensor = matrix2;
    using rank_four_tensor = matrix3;
};

template <>
struct traits<type::solid, bool is_symmetric_ = true>
{
    static auto constexpr value = type::solid;

    /** Solid continuum mechanics has three dimension (x, y, z) */
    static auto constexpr size = 3;

    /** matrix type for the deformation and stress vector */
    using rank_two_tensor = matrix3;

    /** matrix type for a fourth order tensor */
    using rank_four_tensor = std::conditional<is_symmetric_, matrix6, matrix9>;

    using internal_variable_type = internal_variable<rank_two_tensor, rank_four_tensor>;
    using constitutive_type = constitutive_model<internal_variable_type>;

    using mesh_type = fem_mesh<value>;
};
}
