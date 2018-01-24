
#pragma once

#include "AffineMicrosphere.hpp"

namespace neon::mechanical::solid
{
/**
 * NonAffineMicrosphere model computes the Kirchhoff stress and the material
 * tangent for the non-affine microsphere model \cite Miehe2004.  This model includes
 * the interaction of the polymer chain with the neighbouring (forest) chains
 * through the inclusion of the tube model.  This is effectively an extension
 * of the AffineMicrosphere constitutive model.
 */
class NonAffineMicrosphere : public AffineMicrosphere
{
public:
    /**
     * @param variables Reference to internal state variable store
     * @param material_data Json object with material data
     */
    explicit NonAffineMicrosphere(std::shared_ptr<InternalVariables>& variables,
                                  Json::Value const& material_data,
                                  unit_sphere_quadrature::Rule const rule);

    virtual void update_internal_variables(double const time_step_size) override;

protected:
    /**
     * Compute the non-affine stretch by numerical quadrature over the unit sphere
     * \f{align*}{
        \lambda &= \left[ \sum_{i=1}^m (\bar{\lambda}_i)^p w_i \right]^{1/p}
      \f}
     */
    double compute_nonaffine_stretch(matrix3 const& F_deviatoric) const;

    /**
     * Evaluate the h tensor by numerical quadrature over the unit sphere
     * \f{align*}{
        \mathbf{h} &= \sum_{i=1}^m (\bar{\lambda}_i)^{p - 2}
                      \mathbf{t}_i\otimes\mathbf{t}_i w_i
      \f}
     */
    matrix3 compute_h_tensor(matrix3 const& F_deviatoric) const;

    /**
     * Evaluate the H tensor by numerical quadrature over the unit sphere
     * \f{align*}{
        \mathbf{H} &= (p - 2) \sum_{i=1}^m (\bar{\lambda}_i)^{p-4}
                      \mathbf{t}_i\otimes\mathbf{t}_i\otimes\mathbf{t}_i\otimes\mathbf{t}_i w_i
      \f}
     */
    matrix6 compute_H_tensor(matrix3 const& F_deviatoric) const;

    /**
     * Evaluate the k tensor by numerical quadrature over the unit sphere
     * \f{align*}{
        \mathbf{k} &= q \sum_{i=1}^m (\bar{\nu}_i)^{q-2}
                      \mathbf{n}_i\otimes\mathbf{n}_i w_i
      \f}
     * where \f$ q \f$ is the non-affine tube parameter and \f$ \mathbf{n} \f$
     * is the deformed normal.
     * \sa non_affine_tube_parameter
     * @param F_deviatoric Deviatoric part deformation gradient
     */
    matrix3 compute_k_tensor(matrix3 const& F_deviatoric) const;

    /**
     * Evaluate the K tensor by numerical quadrature over the unit sphere
     * \f{align*}{
        \mathbf{K} &= q(q-2)\sum_{i=1}^m (\bar{\nu}_i)^{q-4}
                      \mathbf{n}_i\otimes\mathbf{n}_i\otimes\mathbf{n}_i\otimes\mathbf{n}_iw_i
      \f}
     * where \f$ q \f$ is the non-affine tube parameter.
     * \sa non_affine_tube_parameter
     */
    matrix6 compute_K_tensor(matrix3 const& F_deviatoric) const;

    /**
     * Evaluate the G tensor by numerical quadrature over the unit sphere
     * \f{align*}{
        \mathbf{G} &= 2q \sum_{i=1}^m (\bar{\nu}_i)^{q - 2}
        sym[\mathbf{g}^{-1} \odot \mathbf{n} \otimes \mathbf{n} + \mathbf{n} \otimes
    \mathbf{n} \odot \mathbf{g}^{-1}] w_i
      \f}
     * where \f$ q \f$ is the non-affine tube parameter.
     * \sa non_affine_tube_parameter
     * \sa compute_o_dot_product
     */
    matrix6 compute_G_tensor(matrix3 const& F_deviatoric) const;

    /**
     * Compute the o dot product for the material tangent matrix for the
     * non-affine model using the formula
     * \f{align*}{
         sym[\mathbf{g}^{-1} \odot \mathbf{n} \otimes \mathbf{n} + \mathbf{n} \otimes
     \mathbf{n} \odot \mathbf{g}^{-1}] & \f}
     * where \f$\mathbf{g}\f$ is the metric tensor (identity) and \f$\mathbf{n}\f$
     * is the deformed normal vector
     * \sa compute_G_tensor
     */
    matrix6 compute_o_dot_product(vector3 const& n) const;

    vector3 deformed_normal(matrix3 const& F_unimodular, vector3 const& surface_vector) const
    {
        return F_unimodular.inverse().transpose() * surface_vector;
    }

    auto compute_area_stretch(vector3 const& deformed_normal) const
    {
        return deformed_normal.norm();
    }

private:
    MicromechanicalElastomer material; //!< Material with micromechanical parameters

    double non_affine_stretch_parameter{9.0}; //!< Three-dimensional locking characteristics (p)
    double effective_tube_geometry{9.0};      //!< Additional constraint stiffness (U)
    double non_affine_tube_parameter{1.0};    //!< Shape of constraint stress (q)
};

inline matrix6 NonAffineMicrosphere::compute_o_dot_product(vector3 const& n) const
{
    // clang-format off
    return (matrix6() << 2.0 * n(0) * n(0),               0.0,               0.0,                               0.0,                       n(0) * n(2), n(0) * n(1),       //
                                       0.0, 2.0 * n(1) * n(1),               0.0,                       n(1) * n(2),                               0.0, n(0) * n(1),       //
                                       0.0,               0.0, 2.0 * n(2) * n(2),                       n(1) * n(2),                       n(0) * n(2), 0.0,               //
                                       0.0,       n(1) * n(2),       n(1) * n(2), 0.5 * (n(1) * n(1) + n(2) * n(2)),                 0.5 * n(0) * n(1), 0.5 * n(0) * n(2), //
                               n(0) * n(2),               0.0,       n(0) * n(2),                 0.5 * n(0) * n(1), 0.5 * (n(0) * n(0) + n(2) * n(2)), 0.5 * n(2) * n(1), //
                               n(0) * n(1),       n(0) * n(1),               0.0,                 0.5 * n(0) * n(2),                 0.5 * n(2) * n(1), 0.5 * (n(0) * n(0) + n(1) * n(1))).finished();
    // clang-format on
}
}
