
#pragma once

#include "constitutive/variable_types.hpp"

#include <optional>

namespace neon
{
namespace mechanical
{
std::optional<variable::scalar> is_scalar(std::string const& name);

std::optional<variable::vector> is_vector(std::string const& name);

std::optional<variable::second> is_second_order_tensor(std::string const& name);

//     // clang-format off
//     scalar_map_t const scalar_map{{"AccumulatedPlasticStrain", variable_type::scalar::effective_plastic_strain},
//                                   {"VonMisesStress", variable_type::scalar::von_mises_stress},
//                                   {"Damage", variable_type::scalar::damage},
//                                   {"ActiveShearModulus", variable_type::scalar::active_shear_modulus},
//                                   {"InactiveShearModulus", variable_type::scalar::inactive_shear_modulus},
//                                   {"ActiveSegments", variable_type::scalar::active_segments},
//                                   {"InactiveSegments", variable_type::scalar::inactive_segments},
//                                   {"ReductionFactor", variable_type::scalar::reduction_factor},
//                                   {"EnergyReleaseRate", variable_type::scalar::energy_release_rate}};
//
//     tensor_map_t const tensor_map{{"CauchyStress", variable_type::second::cauchy_stress},
//                                   {"LinearisedStrain", variable_type::second::linearised_strain},
//                                   {"LinearisedPlasticStrain", variable_type::second::linearised_plastic_strain},
//                                   {"DeformationGradient", variable_type::second::deformation_gradient},
//                                   {"DisplacementGradient", variable_type::second::displacement_gradient},
//                                   {"KinematicHardening", variable_type::second::kinematic_hardening},
//                                   {"BackStress", variable_type::second::back_stress}};
//     // clang-format on
}
}
