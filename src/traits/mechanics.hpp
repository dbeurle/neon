
#pragma once

#include <cstdint>

namespace neon::mechanical
{
enum class type : uint8_t { plane_stress, plane_strain, axisymmetric, solid, shell, beam };

template <type T, bool is_symmetric_>
struct traits;

template <>
struct traits<type::plane_strain, bool is_symmetric_ = true>
{
    static auto constexpr value{type::plane_strain};
    static auto constexpr size{2};
    static auto constexpr is_symmetric{is_symmetric_};

    /** Plane stress mechanics has two unknowns (u_x, u_y) */
    static std::array<std::int32_t, 3> constexpr dofs{0, 1};

    using rank_two_tensor = matrix2;
    using rank_four_tensor = matrix3;
};

template <>
struct traits<type::plane_stress, bool is_symmetric_ = true>
{
    static auto constexpr value = type::plane_stress;
    static auto constexpr size{2};
    static auto constexpr is_symmetric = is_symmetric_;

    /** Plane stress mechanics has two unknowns (u_x, u_y) */
    static std::array<std::int32_t, 3> constexpr dofs{0, 1};

    using rank_two_tensor = matrix2;
    using rank_four_tensor = matrix3;
};

template <>
struct traits<type::solid, bool is_symmetric_ = true>
{
    static auto constexpr value{type::solid};
    static auto constexpr size{3};
    static auto constexpr is_symmetric{is_symmetric_};

    /** Solid mechanics has three unknowns (u_x, u_y, u_z) */
    static std::array<std::int32_t, 3> constexpr dofs{0, 1, 2};

    /** matrix type for the deformation and stress vector */
    using rank_two_tensor = matrix3;
    /** matrix type for a fourth order tensor */
    using rank_four_tensor = std::conditional<is_symmetric, matrix6, matrix9>::type;
};
}
