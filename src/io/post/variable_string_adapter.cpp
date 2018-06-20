
#include "variable_string_adapter.hpp"

namespace neon::variable
{
std::optional<variable::scalar> is_scalar(std::string const& name)
{
    if (name == "active_shear_modulus")
    {
        return std::make_optional(scalar::active_shear_modulus);
    }
    else if (name == "inactive_shear_modulus")
    {
        return std::make_optional(scalar::inactive_shear_modulus);
    }
    else if (name == "active_segments")
    {
        return std::make_optional(scalar::active_segments);
    }
    else if (name == "inactive_segments")
    {
        return std::make_optional(scalar::inactive_segments);
    }
    else if (name == "reduction_factor")
    {
        return std::make_optional(scalar::reduction_factor);
    }
    else if (name == "von_mises_stress")
    {
        return std::make_optional(scalar::von_mises_stress);
    }
    else if (name == "effective_plastic_strain")
    {
        return std::make_optional(scalar::effective_plastic_strain);
    }
    else if (name == "DetF0")
    {
        return std::make_optional(scalar::DetF0);
    }
    else if (name == "DetF")
    {
        return std::make_optional(scalar::DetF);
    }
    else if (name == "damage")
    {
        return std::make_optional(scalar::damage);
    }
    else if (name == "energy_release_rate")
    {
        return std::make_optional(scalar::energy_release_rate);
    }
    else if (name == "torsional_stiffness")
    {
        return std::make_optional(scalar::torsional_stiffness);
    }
    else if (name == "axial_stiffness")
    {
        return std::make_optional(scalar::axial_stiffness);
    }
    else if (name == "shear_area_1")
    {
        return std::make_optional(scalar::shear_area_1);
    }
    else if (name == "shear_area_2")
    {
        return std::make_optional(scalar::shear_area_2);
    }
    else if (name == "cross_sectional_area")
    {
        return std::make_optional(scalar::cross_sectional_area);
    }
    else if (name == "second_moment_area_1")
    {
        return std::make_optional(scalar::second_moment_area_1);
    }
    else if (name == "second_moment_area_2")
    {
        return std::make_optional(scalar::second_moment_area_2);
    }
    return std::nullopt;
}

std::optional<variable::vector> is_vector(std::string const& name) { return std::nullopt; }

std::optional<variable::second> is_second_order_tensor(std::string const& name)
{
    if (name == "cauchy_stress")
    {
        return std::make_optional(second::cauchy_stress);
    }
    else if (name == "kirchhoff_stress")
    {
        return std::make_optional(second::kirchhoff_stress);
    }
    else if (name == "linearised_strain")
    {
        return std::make_optional(second::linearised_strain);
    }
    else if (name == "linearised_plastic_strain")
    {
        return std::make_optional(second::linearised_plastic_strain);
    }
    else if (name == "hencky_strain_elastic")
    {
        return std::make_optional(second::hencky_strain_elastic);
    }
    else if (name == "deformation_gradient")
    {
        return std::make_optional(second::deformation_gradient);
    }
    else if (name == "deformation_gradient_plastic")
    {
        return std::make_optional(second::deformation_gradient_plastic);
    }
    else if (name == "displacement_gradient")
    {
        return std::make_optional(second::displacement_gradient);
    }
    else if (name == "green_lagrange")
    {
        return std::make_optional(second::green_lagrange);
    }
    else if (name == "back_stress")
    {
        return std::make_optional(second::back_stress);
    }
    else if (name == "kinematic_hardening")
    {
        return std::make_optional(second::kinematic_hardening);
    }
    else if (name == "conductivity")
    {
        return std::make_optional(second::conductivity);
    }
    else if (name == "bending_stiffness")
    {
        return std::make_optional(second::bending_stiffness);
    }
    else if (name == "shear_stiffness")
    {
        return std::make_optional(second::shear_stiffness);
    }
    return std::nullopt;
}
}
