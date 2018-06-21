
#pragma once

#include <variant>

/// Define namespace for variable enumeration types
namespace neon::variable
{
enum class nodal : short { displacement, rotation, reaction_force };

/// Names for scalar values
enum class scalar : short {
    /// Active chains per unit volume
    active_shear_modulus,
    /// Inactive chains per unit volume
    inactive_shear_modulus,
    /// Active number of segments per chain
    active_segments,
    /// Inactive number of segments per chain
    inactive_segments,
    /// Reduction factor for network decay
    reduction_factor,
    /// Von Mises equivalent stress
    von_mises_stress,
    /// Accumulated (effective) plastic strain
    effective_plastic_strain,
    /// Reference Jacobian determinant
    DetF0,
    /// Updated Jacobian determinant
    DetF,
    /// Scalar damage variable
    damage,
    /// Scalar energy release rate
    energy_release_rate,
    /// Beam torsional stiffness
    torsional_stiffness,
    /// Beam axial stiffness
    axial_stiffness,
    /// Beam section shear area (one)
    shear_area_1,
    /// Beam section shear area (two)
    shear_area_2,
    /// Section cross sectional area
    cross_sectional_area,
    /// Beam second moment of area (one)
    second_moment_area_1,
    /// Beam second moment of area (two)
    second_moment_area_2
};

/// Names for vector values using a std::vector type
enum class vector : short {
    /// Primary ageing shear modulus
    first_ageing_moduli,
    /// Secondary ageing shear modulus
    second_ageing_moduli,
    /// Tertiary ageing shear modulus
    third_ageing_moduli,
    first_previous,
    second_previous,
    third_previous
};

/// Second order tensor internal variables types
enum class second : short {
    /// Cauchy stress
    cauchy_stress,
    /// Kirchhoff stress
    kirchhoff_stress,
    /// First Piola-Kirchhoff stress (PK1)
    piola_kirchhoff1,
    /// Second Piola-Kirchhoff stress (PK2)
    piola_kirchhoff2,
    /// Linearised (infinitesimal) total strain
    linearised_strain,
    /// Linearised (infinitesimal) plastic strain
    linearised_plastic_strain,
    /// Hencky elastic strain
    hencky_strain_elastic,
    /// Deformation gradient (F)
    deformation_gradient,
    /// Plastic deformation gradient (Fp)
    deformation_gradient_plastic,
    /// Displacement gradient (H)
    displacement_gradient,
    /// Green-Lagrange strain (E)
    green_lagrange,
    /// Back stress (for hardening)
    back_stress,
    /// Kinematic hardening
    kinematic_hardening,
    /// Conductivity tensor
    conductivity,
    /// Beam bending stiffness
    bending_stiffness,
    /// Beam shear stiffness
    shear_stiffness
};

/// Fourth order tensor types
enum class fourth : short {
    /// Material tangent operator
    tangent_operator
};

using types = std::variant<scalar, vector, second, fourth, nodal>;
}
