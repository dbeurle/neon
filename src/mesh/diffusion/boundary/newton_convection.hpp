
#pragma once

#include "mesh/generic/neumann.hpp"

#include "interpolations/shape_function.hpp"

namespace neon::diffusion::boundary
{
/**
 * newton_convection is responsible for computing the additional term to the
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
class newton_convection : public surface_load<surface_interpolation>
{
public:
    /**
     * Construct the boundary condition
     * @param sf A surface interpolation
     * @param nodal_connectivity
     * @param times A list of times and values
     * @param heat_flux A list of heat fluxes
     * @param external_temperature A list of heat transfer coefficients
     */
    explicit newton_convection(std::unique_ptr<surface_interpolation>&& sf,
                               indices const& nodal_connectivity,
                               indices const& dof_list,
                               std::shared_ptr<material_coordinates>& mesh_coordinates,
                               json const& times,
                               json const& heat_flux,
                               json const& heat_transfer_coefficient);

    /**
     * Compute the element stiffness matrix contributing to the mixed boundary condition
       \f{align*}{
         k_{ab} &= \int_{\Gamma} N_a \lambda N_b d\Gamma
       \f}
     */
    [[nodiscard]] std::pair<index_view, matrix> external_stiffness(std::int64_t const element,
                                                                   double const load_factor) const;

protected:
    std::vector<std::pair<double, double>> stiffness_time_data;
};
}
