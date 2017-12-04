
#pragma once

#include "mesh/generic/Neumann.hpp"

#include "interpolations/ShapeFunction.hpp"

namespace neon::diffusion
{
/**
 * NewtonConvection is responsible for computing the additional term to the
 * stiffness matrix and right hand side for modelling physics such as the Newton's
 * law of cooling given by
 *
   \f{align*}{
     h &= \lambda u - q_i n_i
   \f}
 * where \f$ h \f$ is the heat flux (\f$ W/(m^2) \f$), \f$ \lambda \f$ is the heat
 * transfer coefficient \f$ W / (m^2 K) \f$ and \f$ q_i n_i \f$ is the heat
 * leaving the surface.  The surrounding temperature \f$T_\infty\f$ is computed by
 * \f$ T_\infty = h / \lambda \f$
 */
class NewtonConvection : public SurfaceLoad<SurfaceInterpolation>
{
public:
    explicit NewtonConvection(std::unique_ptr<SurfaceInterpolation>&& sf,
                              std::vector<List> const& nodal_connectivity,
                              std::vector<List> const& dof_list,
                              std::shared_ptr<MaterialCoordinates>& material_coordinates,
                              Json::Value const& time_history,
                              Json::Value const& heat_transfer_coefficients,
                              Json::Value const& prescribed_heat_fluxes);

    /**
     * Compute the element stiffness matrix contributing to the mixed boundary condition
       \f{align*}{
         k_{ab} &= \int_{\Gamma} N_a \lambda N_b d\Gamma
       \f}
     */
    [[nodiscard]] std::tuple<List const&, Matrix> external_stiffness(int const element,
                                                                     double const load_factor) const;

protected:
    std::vector<std::pair<double, double>> coefficient_time_data;
};
}
