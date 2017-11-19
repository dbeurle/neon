
#pragma once

#include "J2Plasticity.hpp"

#include "numeric/Tensor.hpp"

#include "material/IsotropicElasticPlasticDamage.hpp"

namespace neon::mech::solid
{
class J2PlasticityDamage : public J2Plasticity
{
public:
    J2PlasticityDamage(InternalVariables& variables, Json::Value const& material_data);

    ~J2PlasticityDamage();

    /**
    * Comute the update of the internal variables based on a standard ductile Chaboche damage model
    [Nonlinear computational structural mechanics, Computational methods for plasticity]
    * The evolution equations are the following
              {\dot{\mathbf{\varepsilon}}}^{p} = \dot{\lambda}_{vp} \ \mathbb{\tilde{N}} \
    \frac{1}{1-D} \\
             {\dot{\mathbf{\alpha}}} = \dot{\lambda}_{vp} \ \left( \mathbb{\tilde{N}} -
    \frac{\gamma}{C} \ \beta \right ) \\
             \dot{D} = \dot{\lambda}_{d}\\
             \mathbf{\tilde{N}} = \frac{3}{2} \frac{\mathbf{\tau}}{J_2(\mathbf{\tau})}\\
             \mathbf{\tau}= \frac{\mathbf{\sigma}^{\rm D}}{1-D} - \mathbf{\beta^{\rm D}}\\
    */
    void update_internal_variables(double const time_step_size) override;

    Material const& intrinsic_material() const override { return material; }

    virtual bool is_finite_deformation() const override { return false; }

    virtual bool is_symmetric() const { return false; };

protected:
    Matrix3 compute_cauchy_stress(Matrix3 const& elastic_strain) const;

    /**
     * Performs the radial return algorithm with nonlinear kinematic hardening.
     * This provides the plastic increment required for updating the internal variables
      * This is done by solving the following system of equations
      * \f[ \dot{\lambda}_{vp} - k_{vp} \langle f_{vp} \rangle^{n_{vp}} = 0 \\
     * \dot{\lambda}_{d} - k_{d} \langle f_{d} \rangle^{n_{d}} = 0 \\
     * \mathbf{\sigma}-(1-D)\ \mathbb{C}\ {\mathbf{\varepsilon}}^{e}_{tr} + \mathbb{C}\
     \Delta{\lambda}_{vp} \ \mathbb{\tilde{N}} = \mathbf{0} \\
     * \mathbf{\beta}-C \ \alpha_n - C * \Delta {\lambda}_{vp} \ \left(
     \mathbb{\tilde{N}}-\frac{\gamma}{C} \ \beta \right)  = \mathbf{0} \\
     * D-D_n-\Delta \ {\lambda}_{d} = 0 \\
     * Y - \frac{1}{2} \left( {\mathbf{\varepsilon}}^{e}_{tr} - \Delta {\lambda}_{vp} \
     \mathbb{\tilde{N}}\ \frac{1}{1-D} \right ) : \mathbb{C} : \left(
     {\mathbf{\varepsilon}}^{e}_{tr}
     - \Delta
     {\lambda}_{vp} \ \mathbb{\tilde{N}} \ \frac{1}{1-D} \right ) = 0 \f]
      * Note that all the n+1 indices are ommited for the sake of simplicity
      */
    double perform_radial_return(Matrix3& cauchy_stress,
                                 Matrix3& back_stress,
                                 double& damage_var,
                                 Matrix3& kin_hard,
                                 double& energy_var,
                                 Matrix6& C_list,
                                 double const& delta_t,
                                 Matrix3 const& eps_e_t);
    /**
     * Evaluates the yield function and returns greater than zero if
     * the yield function has been violated
     */
    double evaluate_yield_function(double const von_mises, Matrix3 const& back_stress) const;

    double evaluate_damage_yield_function(double const energy_var) const;

    Matrix3 compute_stress_like_matrix(Matrix6 const& tangent_operator,
                                       Matrix3 const& strain_like) const;
    Vector6 compute_stress_like_vector(Matrix6 const& tangent_operator,
                                       Matrix3 const& strain_like) const;

protected:
    IsotropicElasticPlasticDamage material;
};
} // namespace neon::mech::solid
