
#include "variable_string_adapter.hpp"

#include <stdexcept>

namespace neon::variable
{
types convert(std::string const& name)
{
    if (name == "active_shear_modulus")
    {
        return scalar::active_shear_modulus;
    }
    else if (name == "inactive_shear_modulus")
    {
        return scalar::inactive_shear_modulus;
    }
    else if (name == "active_segments")
    {
        return scalar::active_segments;
    }
    else if (name == "inactive_segments")
    {
        return scalar::inactive_segments;
    }
    else if (name == "reduction_factor")
    {
        return scalar::reduction_factor;
    }
    else if (name == "von_mises_stress")
    {
        return scalar::von_mises_stress;
    }
    else if (name == "effective_plastic_strain")
    {
        return scalar::effective_plastic_strain;
    }
    else if (name == "DetF0")
    {
        return scalar::DetF0;
    }
    else if (name == "DetF")
    {
        return scalar::DetF;
    }
    else if (name == "damage")
    {
        return scalar::damage;
    }
    else if (name == "energy_release_rate")
    {
        return scalar::energy_release_rate;
    }
    else if (name == "torsional_stiffness")
    {
        return scalar::torsional_stiffness;
    }
    else if (name == "axial_stiffness")
    {
        return scalar::axial_stiffness;
    }
    else if (name == "shear_area_1")
    {
        return scalar::shear_area_1;
    }
    else if (name == "shear_area_2")
    {
        return scalar::shear_area_2;
    }
    else if (name == "cross_sectional_area")
    {
        return scalar::cross_sectional_area;
    }
    else if (name == "second_moment_area_1")
    {
        return scalar::second_moment_area_1;
    }
    else if (name == "second_moment_area_2")
    {
        return scalar::second_moment_area_2;
    }
    else if (name == "cauchy_stress")
    {
        return second::cauchy_stress;
    }
    else if (name == "kirchhoff_stress")
    {
        return second::kirchhoff_stress;
    }
    else if (name == "linearised_strain")
    {
        return second::linearised_strain;
    }
    else if (name == "linearised_plastic_strain")
    {
        return second::linearised_plastic_strain;
    }
    else if (name == "hencky_strain_elastic")
    {
        return second::hencky_strain_elastic;
    }
    else if (name == "deformation_gradient")
    {
        return second::deformation_gradient;
    }
    else if (name == "deformation_gradient_plastic")
    {
        return second::deformation_gradient_plastic;
    }
    else if (name == "displacement_gradient")
    {
        return second::displacement_gradient;
    }
    else if (name == "green_lagrange")
    {
        return second::green_lagrange;
    }
    else if (name == "back_stress")
    {
        return second::back_stress;
    }
    else if (name == "kinematic_hardening")
    {
        return second::kinematic_hardening;
    }
    else if (name == "conductivity")
    {
        return second::conductivity;
    }
    else if (name == "bending_stiffness")
    {
        return second::bending_stiffness;
    }
    else if (name == "shear_stiffness")
    {
        return second::shear_stiffness;
    }
    else if (name == "displacement")
    {
        return nodal::displacement;
    }
    else if (name == "rotation")
    {
        return nodal::rotation;
    }
    else if (name == "reaction_force")
    {
        return nodal::reaction_force;
    }
    throw std::domain_error("Name " + name + " is not a valid output variable name\n");
    // dummy return
    return second::cauchy_stress;
}

char const* convert(variable::scalar value)
{
    switch (value)
    {
        case variable::scalar::active_shear_modulus:
            return "active_shear_modulus";
            break;
        case variable::scalar::inactive_shear_modulus:
            return "inactive_shear_modulus";
            break;
        case variable::scalar::active_segments:
            return "active_segments";
            break;
        case variable::scalar::inactive_segments:
            return "inactive_segments";
            break;
        case variable::scalar::reduction_factor:
            return "reduction_factor";
            break;
        case variable::scalar::von_mises_stress:
            return "von_mises_stress";
            break;
        case variable::scalar::effective_plastic_strain:
            return "effective_plastic_strain";
            break;
        case variable::scalar::DetF0:
            return "DetF0";
            break;
        case variable::scalar::DetF:
            return "DetF";
            break;
        case variable::scalar::damage:
            return "damage";
            break;
        case variable::scalar::energy_release_rate:
            return "energy_release_rate";
            break;
        case variable::scalar::torsional_stiffness:
            return "torsional_stiffness";
            break;
        case variable::scalar::axial_stiffness:
            return "axial_stiffness";
            break;
        case variable::scalar::shear_area_1:
            return "shear_area_1";
            break;
        case variable::scalar::shear_area_2:
            return "shear_area_2";
            break;
        case variable::scalar::cross_sectional_area:
            return "cross_sectional_area";
            break;
        case variable::scalar::second_moment_area_1:
            return "second_moment_area_1";
            break;
        case variable::scalar::second_moment_area_2:
            return "second_moment_area_2";
            break;
    }
    return "\0";
}

char const* convert(variable::second value)
{
    switch (value)
    {
        case variable::second::cauchy_stress:
            return "cauchy_stress";
            break;
        case variable::second::kirchhoff_stress:
            return "kirchhoff_stress";
            break;
        case variable::second::piola_kirchhoff1:
            return "piola_kirchhoff1";
            break;
        case variable::second::piola_kirchhoff2:
            return "piola_kirchhoff2";
            break;
        case variable::second::linearised_strain:
            return "linearised_strain";
            break;
        case variable::second::linearised_plastic_strain:
            return "linearised_plastic_strain";
            break;
        case variable::second::hencky_strain_elastic:
            return "hencky_strain_elastic";
            break;
        case variable::second::deformation_gradient:
            return "deformation_gradient";
            break;
        case variable::second::deformation_gradient_plastic:
            return "deformation_gradient_plastic";
            break;
        case variable::second::displacement_gradient:
            return "displacement_gradient";
            break;
        case variable::second::green_lagrange:
            return "green_lagrange";
            break;
        case variable::second::back_stress:
            return "back_stress";
            break;
        case variable::second::kinematic_hardening:
            return "kinematic_hardening";
            break;
        case variable::second::conductivity:
            return "conductivity";
            break;
        case variable::second::bending_stiffness:
            return "bending_stiffness";
            break;
        case variable::second::shear_stiffness:
            return "shear_stiffness";
            break;
    }
    return "\0";
}
}
