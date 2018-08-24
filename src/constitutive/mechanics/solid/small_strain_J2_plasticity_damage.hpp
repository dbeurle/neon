
#pragma once

#include "small_strain_J2_plasticity.hpp"

#include "numeric/tensor_operations.hpp"

#include "material/isotropic_elastic_plastic_damage.hpp"

namespace neon::mechanics::solid
{
/// small_strain_J2_plasticity_damage is responsible for computing the
/// constitutive model of a standard ductile Chaboche damage model without
/// microdefects closure effects. Theoritical details can be found in
/// \cite Neto2011 (Lemaitre’s elastoplastic damage theory),
/// \cite Ladeveze1999 (chapter 2) and
/// \cite Lemaitre2005 (Isotropic Damage without Microdefects Closure).
class small_strain_J2_plasticity_damage : public small_strain_J2_plasticity
{
public:
    small_strain_J2_plasticity_damage(std::shared_ptr<internal_variables_t>& variables,
                                      json const& material_data);

    ~small_strain_J2_plasticity_damage();

    /**
    * Compute the update of the internal variables based on the following evolution equations.
    * \f[          {\dot{\mathbf{\varepsilon}}}^{p} = \dot{\lambda}_{vp} \ \mathbf{\tilde{N}} \
    \frac{1}{1-D} \\
             {\dot{\mathbf{\alpha}}} = \dot{\lambda}_{vp} \ \left( \mathbf{\tilde{N}} -
    \frac{\gamma}{C} \ \beta \right ) \\
             \dot{D} = \dot{\lambda}_{d} \f]
             with \f$ \mathbf{\tilde{N}} = \frac{3}{2} \frac{\mathbf{\tau}}{J_2(\mathbf{\tau})}\f$
    and \f$ \mathbf{\tau}= \frac{\mathbf{\sigma}^{\rm D}}{1-D} - \mathbf{\beta^{\rm D}} \f$
    */
    void update_internal_variables(double const time_step_size) override;

    material_property const& intrinsic_material() const override { return material; }

    virtual bool is_finite_deformation() const override { return false; }

    virtual bool is_symmetric() const override { return false; };

protected:
    /**
     * Performs the radial return algorithm with nonlinear kinematic hardening and damage.
     * This provides the plastic increment required for updating the plastic strain. cauchy_stress,
     back_stress, damage_var, kin_hard, energy_var and C_algorithmic, are updated within this
     routine.
      * This is done by solving the following nonlinear system of equations:
      * \f[ \dot{\lambda}_{vp} - k_{vp} \langle f_{vp} \rangle^{n_{vp}} = 0 \\
     * \dot{\lambda}_{d} - k_{d} \langle f_{d} \rangle^{n_{d}} = 0 \\
     * \mathbf{\sigma}-(1-D)\ \mathbb{C}\ {\mathbf{\varepsilon}}^{e}_{tr} + \mathbb{C}\
     \Delta{\lambda}_{vp} \ \mathbf{\tilde{N}} = \mathbf{0} \\
     * \mathbf{\beta}-C \ \alpha_n - C \ \Delta {\lambda}_{vp} \ \left(
     \mathbf{\tilde{N}}-\frac{\gamma}{C} \ \beta \right)  = \mathbf{0} \\
     * D-D_n-\Delta \ {\lambda}_{d} = 0 \\
     * Y - \frac{1}{2} \left( {\mathbf{\varepsilon}}^{e}_{tr} - \Delta {\lambda}_{vp} \
     \mathbf{\tilde{N}}\ \frac{1}{1-D} \right ) : \mathbb{C} : \left(
     {\mathbf{\varepsilon}}^{e}_{tr}
     - \Delta
     {\lambda}_{vp} \ \mathbf{\tilde{N}} \ \frac{1}{1-D} \right ) = 0 \f]
      * Note that all the n+1 indices are ommited for the sake of simplicity.
      */
    [[nodiscard]] double perform_radial_return(matrix3& cauchy_stress,
                                               matrix3& back_stress,
                                               double& damage_var,
                                               matrix3& kin_hard,
                                               double& energy_var,
                                               matrix6& C_list,
                                               double const& delta_t,
                                               matrix3 const& eps_e_t);

    /// Evaluates the yield function and returns greater than zero if
    /// the yield function has been violated
    [[nodiscard]] double evaluate_yield_function(double const von_mises,
                                                 matrix3 const& back_stress) const;

    [[nodiscard]] double evaluate_damage_yield_function(double const energy_var) const;

    [[nodiscard]] matrix3 compute_stress_like_matrix(matrix6 const& tangent_operator,
                                                     matrix3 const& strain_like) const;

    [[nodiscard]] vector6 compute_stress_like_vector(matrix6 const& tangent_operator,
                                                     matrix3 const& strain_like) const;

protected:
    isotropic_elastic_plastic_damage material;
};
}
